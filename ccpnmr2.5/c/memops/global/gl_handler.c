
/*
======================COPYRIGHT/LICENSE START==========================

gl_handler.c: Part of the CcpNmr Analysis program

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
#include "gl_handler.h"

#ifndef IGNORE_GL

#ifdef WIN64
#include <windows.h>
#include <tkPlatDecls.h>
#else
#include <GL/glx.h>
#endif /* end WIN32 */

#include <GL/gl.h>
#ifdef WIN64
#include <GL/glu.h>
#endif

#ifdef GLUT_IN_OWN_DIR
#include "glut.h"
#else
#include "GL/glut.h"
#endif /* end GLUT_IN_OWN_DIR */

#include "clipping.h"
#include "utility.h"


//#include <mach-o/dyld.h>
//
//void * MyNSGLGetProcAddress (const char *name)
//{
//    NSSymbol symbol;
//    char *symbolName;
//    symbolName = malloc (strlen (name) + 2); // 1
//    strcpy(symbolName + 1, name); // 2
//    symbolName[0] = '_'; // 3
//    symbol = NULL;
//    if (NSIsSymbolNameDefined (symbolName)) // 4
//        symbol = NSLookupAndBindSymbol (symbolName);
//    free (symbolName); // 5
//    return symbol ? NSAddressOfSymbol (symbol) : NULL; // 6
//
//}

#define  NCOLORS  3

typedef struct _GlutFont
{
    char *name;
    int size;
    void *font;
    Bool isStrokeFont;
}   GlutFont;

static GlutFont glutBitmapFonts[] =
{
#ifdef GLUT_BITMAP_HELVETICA_10
    { "Helvetica", 10, (void *) GLUT_BITMAP_HELVETICA_10, CCPN_FALSE },
#endif
#ifdef GLUT_BITMAP_HELVETICA_12
    { "Helvetica", 12, (void *) GLUT_BITMAP_HELVETICA_12, CCPN_FALSE },
#endif
#ifdef GLUT_BITMAP_HELVETICA_18
    { "Helvetica", 18, (void *) GLUT_BITMAP_HELVETICA_18, CCPN_FALSE },
#endif

#ifdef GLUT_BITMAP_TIMES_ROMAN_10
    { "Times", 10, (void *) GLUT_BITMAP_TIMES_ROMAN_10, CCPN_FALSE },
#endif
#ifdef GLUT_BITMAP_TIMES_ROMAN_24
    { "Times", 24, (void *) GLUT_BITMAP_TIMES_ROMAN_24, CCPN_FALSE },
#endif

#ifdef GLUT_BITMAP_8_BY_13
    { "Courier", 13, (void *) GLUT_BITMAP_8_BY_13, CCPN_FALSE },
#endif
#ifdef GLUT_BITMAP_9_BY_15
    { "Courier", 15, (void *) GLUT_BITMAP_9_BY_15, CCPN_FALSE },
#endif
};

static int nbitmapFonts = sizeof(glutBitmapFonts) / sizeof(glutBitmapFonts[0]);

static GlutFont glutStrokeFonts[] =
{
#ifdef GLUT_STROKE_ROMAN
    { "Roman", 1, (void *) GLUT_STROKE_ROMAN, CCPN_TRUE },
#endif
#ifdef GLUT_STROKE_ROMAN
    { "MonoRoman", 1, (void *) GLUT_STROKE_MONO_ROMAN, CCPN_TRUE },
#endif
};

static int nstrokeFonts = sizeof(glutStrokeFonts) / sizeof(glutStrokeFonts[0]);

typedef struct Gl_handler_p
{
    Tk_Window tk_win;
    Tcl_Interp *interp;
    Display *display;
#ifdef WIN64
    HWND hWnd;
    HGLRC hGLRC;
#else
    XVisualInfo *visual;
    GLXContext context;
#endif
/*
    GC gc;
*/
    int width;
    int height;
    float background[NCOLORS];
    float x0, y0, x1, y1;
    float sx, sy, tx, ty; /* scalings: not sure need these */
    GLuint main_call_list; /* for main drawing */
    GLuint xor_call_list; /* for XOR drawing */
    int line_style;
    GlutFont *glutFont;
    Bool is_double_buffer;
}   *Gl_handler_p;

/* looks like cannot use X if you are using GL double buffers
   because X draws into the screen buffer rather than the GL buffer */
/*
#define  CONVERT_X(x)  ((int) (gl_handler->width * (gl_handler->sx*(x) + gl_handler->tx)))
#define  CONVERT_Y(y)  ((int) (gl_handler->height * (gl_handler->sy*(y) + gl_handler->ty)))
#define  CONVERT_X(x)  (gl_handler->sx*(x) + gl_handler->tx)
#define  CONVERT_Y(y)  (gl_handler->sy*(y) + gl_handler->ty)
*/

#define  DEFAULT_LINEWIDTH  1.0

static GLushort determinePattern(int dash_length, int gap_length)
{
    int i, nbits = 16, d = dash_length+gap_length;
    GLushort basic, pattern;

    if (gap_length == 0)
    {
	pattern = (1 << nbits) - 1;
    }
    else if (dash_length == 0)
    {
	pattern = 0;
    }
    else
    {
	pattern = 0;
	basic = (1 << dash_length) - 1;;
	for (i = 0; i < nbits; i += d)
	    pattern |= basic << i;

	pattern &= (1 << nbits) - 1;
    }

/*
    printf("determinePattern: dash_length=%d, gap_length=%d, pattern=%d %o\n",
	dash_length, gap_length, pattern, pattern);
*/

    return pattern;
}

static void gl_start_draw(Generic_ptr data)
{
    Gl_handler gl_handler = (Gl_handler) data;

    make_current_gl_handler(gl_handler);
}

static void gl_end_draw(Generic_ptr data)
{
/* TODO */
}

static void gl_new_draw_range(Generic_ptr data, float x0, float y0,
					float x1, float y1, Bool clip)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) data;

    gl_handler_p->x0 = MIN(x0, x1);
    gl_handler_p->y0 = MIN(y0, y1);
    gl_handler_p->x1 = MAX(x0, x1);
    gl_handler_p->y1 = MAX(y0, y1);
}

static void gl_draw_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) data;
    static int dash_length = 2, gap_length = 2;

    if ((gl_handler_p != NULL) && (gl_handler_p->line_style == DASHED_LINE_STYLE))
	draw_dash_line_gl_handler((Gl_handler) data, x0, y0, x1, y1,
					dash_length, gap_length);
    else
	draw_line_gl_handler((Gl_handler) data, x0, y0, x1, y1);
}

static void gl_draw_clipped_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) data;
    static int dash_length = 2, gap_length = 2;

    if ((gl_handler_p != NULL) && (gl_handler_p->line_style == DASHED_LINE_STYLE))
	draw_clipped_dash_line_gl_handler((Gl_handler) data, x0, y0, x1, y1,
					dash_length, gap_length);
    else
	draw_clipped_line_gl_handler((Gl_handler) data, x0, y0, x1, y1);
}

static void gl_draw_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_polyline_gl_handler((Gl_handler) data, polyline);
}

static void gl_draw_clipped_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_clipped_polyline_gl_handler((Gl_handler) data, polyline);
}

static void gl_draw_text(Generic_ptr data, CcpnString text, float x, float y,
							float a, float b)
{
    draw_text_gl_handler((Gl_handler) data, text, x, y, a, b);
}

static void gl_set_draw_color(Generic_ptr data, float *color)
{
    set_color_gl_handler((Gl_handler) data, color);
}

static void gl_set_draw_font(Generic_ptr data, CcpnString name, int size)
{
    set_font_gl_handler((Gl_handler) data, name, size);
}

static void gl_set_line_style(Generic_ptr data, int line_style)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) data;

    gl_handler_p->line_style = line_style;
}

static void gl_set_line_width(Generic_ptr data, float line_width)
{
    set_line_width_gl_handler((Gl_handler) data, line_width);
}

static void gl_fill_circle(Generic_ptr data, float x, float y, float r)
{
    fill_circle_gl_handler((Gl_handler) data, x, y, r);
}

static void gl_fill_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    fill_ellipse_gl_handler((Gl_handler) data, x, y, rx, ry);
}

static void gl_draw_circle(Generic_ptr data, float x, float y, float r)
{
    draw_circle_gl_handler((Gl_handler) data, x, y, r);
}

static void gl_draw_ellipse(Generic_ptr data, float x, float y, float rx, float ry)
{
    draw_ellipse_gl_handler((Gl_handler) data, x, y, rx, ry);
}

static void gl_get_background(Generic_ptr data, float *color)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) data;

    COPY_VECTOR(color, gl_handler_p->background, NCOLORS);
}

static void gl_get_region(Generic_ptr data, float *x0, float *y0,
                                                    float *x1, float *y1)
{
    get_region_gl_handler((Gl_handler) data, x0, y0, x1, y1);
}

static void gl_get_text_size(Generic_ptr data, CcpnString text, float *w, float*h)
{
    get_text_size_gl_handler((Gl_handler) data, text, w, h);
}

static void gl_fill_triangle(Generic_ptr data, float x0, float y0,
		                                float x1, float y1, float x2, float y2)
{
    fill_triangle_gl_handler((Gl_handler) data, x0, y0, x1, y1, x2, y2);
}

static Drawing_funcs drawing_funcs = { SCREEN_DRAWING,
	gl_start_draw, gl_end_draw, gl_new_draw_range,
	gl_draw_line, gl_draw_clipped_line,
	gl_draw_polyline, gl_draw_clipped_polyline, gl_draw_text,
	gl_set_draw_color, gl_set_draw_font, gl_set_line_style,
        gl_set_line_width, gl_fill_circle, gl_fill_ellipse,
	gl_draw_circle, gl_draw_ellipse, gl_get_background,
	gl_get_region, gl_get_text_size, gl_fill_triangle };

Drawing_funcs *gl_drawing_funcs(void)
{
    return &drawing_funcs;
}

static void create_call_lists(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    //GLuint call_list = glGenLists(1);
    gl_handler_p->main_call_list = glGenLists(1);
    gl_handler_p->xor_call_list = glGenLists(1);

/* not sure if below are needed */
/*
    glNewList(call_list, GL_COMPILE);
    glEndList();
*/

/*
    printf("create_call_lists: %x %x\n", gl_handler_p->main_call_list, gl_handler_p->xor_call_list);
*/
}

static void delete_call_lists(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    //GLuint call_list = gl_handler_p->call_list;
/*
    printf("delete_call_list: %x\n", call_list);
*/
    glDeleteLists(gl_handler_p->main_call_list, 1);
    glDeleteLists(gl_handler_p->xor_call_list, 1);
}

static Bool init_gl_handler(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    int i;

    if (!make_current_gl_handler(gl_handler))
      return CCPN_FALSE;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

/*
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
*/

    glEnable(GL_COLOR_LOGIC_OP);
/*
    glLogicOp(GL_EQUIV);
*/
    glLogicOp(GL_COPY);

    for (i = 0; i < NCOLORS; i++)
	gl_handler_p->background[i] = 1;  /* default = white */

/*
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
*/

    glPushMatrix();

    /* arbitrary */
    gl_handler_p->width = 1;
    gl_handler_p->height = 1;
    gl_handler_p->sx = 1.0;
    gl_handler_p->sy = 1.0;
    gl_handler_p->tx = 0.0;
    gl_handler_p->ty = 0.0;
    gl_handler_p->line_style = NORMAL_LINE_STYLE;

    create_call_lists(gl_handler);

    return CCPN_TRUE;
}

static void end_gl_handler(Gl_handler gl_handler)
{
/*  18 Oct 2004: below causes crashes on some Linux boxes
    and cannot see why really need to do this.
*/
    if (!make_current_gl_handler(gl_handler))
      return;

/*
    glPopMatrix();

    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
*/
/*  But assume ought to delete call list... */
    delete_call_lists(gl_handler);
}

/*
 *  Contents of win32func.c where here - dan
 */

Gl_handler new_gl_handler(Tcl_Interp *interp, Tk_Window tk_win, Bool direct, Bool dbl_buff )
{
    static Bool first_pass = CCPN_TRUE;
    Gl_handler gl_handler;
    Gl_handler_p gl_handler_p;
    Display *display;
#ifdef WIN64
    HWND hWnd;
    HGLRC hGLRC;
    int  iPixelFormat;
    PIXELFORMATDESCRIPTOR pfd;
/* After painting with a common DC, the ReleaseDC function must be called to release the DC. 
Class and private DCs do not have to be released. ReleaseDC must be called from the same thread 
that called GetDC. The number of DCs is limited only by available memory.
*/	
    HDC hDC; /* not to be set in gl_handler */
    const char * strExt;
#else
    int dummy;
    XVisualInfo *visual;
    GLXContext context;
    //static int dblBuf[] = {GLX_RGBA, GLX_DOUBLEBUFFER, None};
    static int dblBuf[] = {GLX_RGBA, None};
#endif
/*
    static int dblBuf[] = {GLX_RGBA, GLX_DEPTH_SIZE, 8, GLX_DOUBLEBUFFER, None};
*/
/*
    GC gc;
*/
    XGCValues gcv;
#ifdef NEED_GLUT_INIT
    int argc = 0;
    char *argv = NULL;
#endif

    Tk_MakeWindowExist(tk_win);
    display = Tk_Display(tk_win);

    /* below seems to crash GL on some systems (at least OSX...) */
    /* glGetError(); */ /* remove existing errors */

#ifdef WIN64

    memset( &pfd, 0, sizeof(PIXELFORMATDESCRIPTOR) );
    pfd.nSize       = sizeof(PIXELFORMATDESCRIPTOR);   // size of this pfd 
    pfd.nVersion    = 1;                               // version number 
    pfd.dwFlags     = PFD_DRAW_TO_WINDOW |             // support window 
                      PFD_SUPPORT_OPENGL |             // support OpenGL
                      PFD_GENERIC_ACCELERATED |
                      (dbl_buff ? PFD_DOUBLEBUFFER : PFD_DOUBLEBUFFER_DONTCARE );
                                                       // double buffered                                                      
    pfd.iPixelType  = PFD_TYPE_RGBA;                   // RGBA type 
    pfd.cColorBits  = 24;                              // 24-bit color depth 
    pfd.cDepthBits  = 32;                              // 32-bit z-buffer 
    pfd.iLayerType  = PFD_MAIN_PLANE;                  // main layer 

    hWnd = Tk_GetHWND(Tk_WindowId(tk_win));

    /*important to get and release hDC*/
    if((hDC = GetDC(hWnd)) == NULL)
        printf("new_gl_handler: GetDC failed\n");

    /* get the best available match of pixel format for the device context */
    if( (iPixelFormat = ChoosePixelFormat(hDC, &pfd)) == 0) 
        printf("new_gl_handler: ChoosePixelFormat failed\n"); 

    /* make that the pixel format of the device context */
    if( SetPixelFormat(hDC, iPixelFormat, &pfd) == FALSE)
        printf("new_gl_handler: SetPixelFormat failed: %d\n", GetLastError());

    if( (hGLRC = wglCreateContext(hDC)) == NULL)
        printf("new_gl_handler: wglCreateContext failed\n"); 

    if( wglMakeCurrent(hDC, hGLRC) == FALSE)
        printf("new_gl_handler: wglMakeCurrent failed\n");

    /*release hDC*/
    if(ReleaseDC(hWnd, hDC) == 0) printf("new_gl_handler: ReleaseDC failed\n");

    /* don't use Tk_GetGC if gc to be dynamically altered... */
/*
    gc = Tk_GetGC(tk_win, GCForeground | GCBackground, &gcv);
    gc = XCreateGC(display, Tk_WindowId(tk_win), GCForeground | GCBackground, &gcv);
*/	
    if (!(glGetError() == GL_NO_ERROR))
      return NULL;

    /* no need for anti-aliasing */
    glDisable( GL_LINE_SMOOTH );

#else

//    //~~~~~~~~~~~~~~~~~
//    // NOTE:ED - code to test the glX binding
//    printf("glXQueryVersion\n");
//    typedef bool (*glXQueryVersionProcPtr)(Display * dpy, int * Major, int * Minor);
//    glXQueryVersionProcPtr pfglXQueryVersion = NULL;
//    pfglXQueryVersion = (glXQueryVersionProcPtr) MyNSGLGetProcAddress("glXQueryVersion");
//    GLint majorGLX, minorGLX = 0;
//    if (!glXQueryVersion(display, &majorGLX, &minorGLX))
//        if (majorGLX <= 1 && minorGLX < 2) {
//            printf("GLX 1.2 or greater is required.\n");
//            return NULL;
//        }
//    printf("GLX version %i.%i\n", majorGLX, minorGLX);
//    //~~~~~~~~~~~~~~~~~

    if (!glXQueryExtension(display, &dummy, &dummy))
    {
      printf("glXQueryExtension failed\n");
      printf("If Linux, check that glx (and dri?) enabled in /etc/X11 config file\n");
      return NULL; /* check that OpenGL's GLX extension supported */
    }

    visual = glXChooseVisual(display, DefaultScreen(display), dblBuf);
    if (!visual)
    {
      printf("glXChooseVisual failed\n");
      printf("This implies that your glx does not support double buffering.");
      return NULL;
    }

    context = glXCreateContext(display, visual, None, direct ? GL_TRUE : GL_FALSE);

/*
    if (!context || !(glGetError() == GL_NO_ERROR))
*/
    if (!context)
    {
      printf("glXCreateContext failed\n");
      printf("If environment.txt has GL_FLAG = -DUSE_GL_TRUE, switch to ..._FALSE, and vice versa.  Then in the environment.txt directory do 'make clean' and 'make'\n");
      return NULL;
    }

    gcv.foreground = BlackPixelOfScreen(DefaultScreenOfDisplay(display));
    gcv.background = WhitePixelOfScreen(DefaultScreenOfDisplay(display));
/*
    gc = XCreateGC(display, Tk_WindowId(tk_win), GCForeground | GCBackground, &gcv);
*/
#endif

    MALLOC_NEW(gl_handler_p, struct Gl_handler_p, 1);

    gl_handler_p->tk_win = tk_win;
    gl_handler_p->interp = interp;
    gl_handler_p->display = display;
#ifdef WIN64
    gl_handler_p->hWnd = hWnd;
    gl_handler_p->hGLRC = hGLRC;
#else
    gl_handler_p->visual = visual;
    gl_handler_p->context = context;
#endif
/*
    gl_handler_p->gc = gc;
*/
    if (nbitmapFonts > 0)
    	gl_handler_p->glutFont = glutBitmapFonts;
    else if (nstrokeFonts > 0)
    	gl_handler_p->glutFont = glutStrokeFonts;
    else
    	gl_handler_p->glutFont = NULL;
    gl_handler_p->is_double_buffer = CCPN_TRUE;

    gl_handler = (Gl_handler) gl_handler_p;

    if (!init_gl_handler(gl_handler))
    {
      FREE(gl_handler_p, struct Gl_handler_p);

      return NULL;
    }

    if (first_pass)
    {
    	first_pass = CCPN_FALSE;
	    printf("GL_VENDOR = %s\n", glGetString(GL_VENDOR));
	    printf("GL_RENDERER = %s\n", glGetString(GL_RENDERER));
	    printf("GL_VERSION = %s\n", glGetString(GL_VERSION));
/*
	printf("GL_EXTENSIONS = %s\n", glGetString(GL_EXTENSIONS));
*/
#ifdef NEED_GLUT_INIT
        glutInit(&argc, &argv);
#endif
    }

/*
    printf("new_gl_handler %x\n", gl_handler);
*/
    return gl_handler;
}

void delete_gl_handler(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

/*
    printf("delete_gl_handler %x\n", gl_handler);
*/

    if (gl_handler)
    {
        end_gl_handler(gl_handler);

#ifdef WIN64
        if (gl_handler_p->hGLRC)
        {
            if( wglMakeCurrent(NULL, NULL) == FALSE)
                printf("delete_gl_handler: wglMakeCurrent failed\n");
            if( wglDeleteContext(gl_handler_p->hGLRC) == FALSE)
                printf("delete_gl_handler: wglDeleteContext failed\n");
        }
#else
 	glXDestroyContext(gl_handler_p->display, gl_handler_p->context);
	if (gl_handler_p->visual)
	    XFree(gl_handler_p->visual);
#endif
/* 	
	if (gl_handler_p->gc)
	    XFreeGC(gl_handler_p->display, gl_handler_p->gc);
*/
	FREE(gl_handler_p, struct Gl_handler_p);
    }
}

void resize_gl_handler(Gl_handler gl_handler, int width, int height)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
/*  does below mean mapped onto (0, 1)???  */
/*
    GLdouble left = 0, right = width, bottom = 0, top = height;
*/
    GLdouble left = 0, right = 1, bottom = 0, top = 1;
    GLdouble nearr = -1, farr = 1; /* extra r required by Windows to avoid key words */

/*
    printf("resize_gl_handler %x\n", gl_handler);
*/

/*
    printf("resize_gl_handler: w = %d, h = %d\n", width, height);
*/
    if (!make_current_gl_handler(gl_handler))
      return;

/*
    XResizeWindow(gl_handler->display, Tk_WindowId(gl_handler->tk_win),
							width, height);
*/

    glPopMatrix();
    glPushMatrix();

    glViewport(0, 0, (GLint) width, (GLint) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(left, right, bottom, top, nearr, farr);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gl_handler_p->width = width;
    gl_handler_p->height = height;
}

void expose_gl_handler(Gl_handler gl_handler, int x, int y, int w, int h)
{
/*
    static GLdouble eqn0[] = { 1.0, 0.0, 0.0, 0.0 };
    static GLdouble eqn1[] = { -1.0, 0.0, 0.0, 0.0 };
    static GLdouble eqn2[] = { 0.0, 1.0, 0.0, 0.0 };
    static GLdouble eqn3[] = { 0.0, -1.0, 0.0, 0.0 };
*/
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLclampf r = gl_handler_p->background[0];
    GLclampf g = gl_handler_p->background[1];
    GLclampf b = gl_handler_p->background[2];

/*
    printf("expose_gl_handler: x = %d, y = %d, w = %d, h = %d\n", x, y, w, h);
    printf("expose_gl_handler: r = %2.1f, g = %2.1f, b = %2.1f\n", r, g, b);
*/
/*
    if (!make_current_gl_handler(gl_handler))
      return;
*/

    glClearColor(r, g, b, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

/*
    eqn0[3] = - (double) x / (double) gl_handler->width;
    eqn1[3] = (double) (x+w) / (double) gl_handler->width;
    eqn2[3] = - (double) y / (double) gl_handler->height;
    eqn3[3] = (double) (y+h) / (double) gl_handler->height;
*/

/*
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
*/

/*
    glClipPlane(GL_CLIP_PLANE0, eqn0);
    glClipPlane(GL_CLIP_PLANE1, eqn1);
    glClipPlane(GL_CLIP_PLANE2, eqn2);
    glClipPlane(GL_CLIP_PLANE3, eqn3);
*/
}

void flush_gl_handler(Gl_handler gl_handler)
{
    glFlush();
}

void map_ranges_gl_handler(Gl_handler gl_handler,
		float x0, float y0, float x1, float y1,
		float a0, float b0, float a1, float b1)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
/*
    static GLdouble eqn0[] = { 1.0, 0.0, 0.0, 0.0 };
    static GLdouble eqn1[] = { -1.0, 0.0, 0.0, 0.0 };
    static GLdouble eqn2[] = { 0.0, 1.0, 0.0, 0.0 };
    static GLdouble eqn3[] = { 0.0, -1.0, 0.0, 0.0 };
*/
    GLdouble sx, sy, tx, ty;

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
    if (!make_current_gl_handler(gl_handler))
      return;
*/

    glPopMatrix();
    glPushMatrix();

/*
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
*/

/*
    eqn0[3] = - x0;
    eqn1[3] = x1;
    eqn2[3] = - y0;
    eqn3[3] = y1;

printf("eqn0 = %6.3f\n", eqn0[3]);
printf("eqn1 = %6.3f\n", eqn1[3]);
printf("eqn2 = %6.3f\n", eqn2[3]);
printf("eqn3 = %6.3f\n", eqn3[3]);

    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);

    glClipPlane(GL_CLIP_PLANE0, eqn0);
    glClipPlane(GL_CLIP_PLANE1, eqn1);
    glClipPlane(GL_CLIP_PLANE2, eqn2);
    glClipPlane(GL_CLIP_PLANE3, eqn3);
*/

    sx = gl_handler_p->sx = (x1 - x0) / (a1 - a0);
    sy = gl_handler_p->sy = (y1 - y0) / (b1 - b0);
    tx = gl_handler_p->tx = (x0*a1-x1*a0) / (a1 - a0);
    ty = gl_handler_p->ty = (y0*b1-y1*b0) / (b1 - b0);

    a0 -= x0 / sx;
    a1 += (1-x1) / sx;
    b0 -= y0 / sy;
    b1 += (1-y1) / sy;
    gl_new_draw_range(gl_handler, a0, b0, a1, b1, CCPN_TRUE);

/*
printf("a0 = %6.3f, b0 = %6.3f\n", a0, b0);
printf("a1 = %6.3f, b1 = %6.3f\n", a1, b1);
printf("x0 = %6.3f, y0 = %6.3f\n", x0, y0);
printf("x1 = %6.3f, y1 = %6.3f\n", x1, y1);
printf("sx = %7.4f, sy = %7.4f\n", sx, sy);
printf("tx = %6.3f, ty = %6.3f\n", tx, ty);
*/

    glTranslated(tx, ty, 0.0);
    glScaled(sx, sy, 1.0);
}

Bool make_current_gl_handler(Gl_handler gl_handler)
{
    Bool ok = 1;
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    /* below seems to crash GL on some systems (at least OSX...) */
    /* glGetError(); */ /* remove existing errors */

#ifdef WIN64
    HDC hDC;
    if((hDC = GetDC(gl_handler_p->hWnd)) == NULL) 
    	printf("make_current_gl_handler: GetDC failed\n");
    if((wglMakeCurrent(hDC, gl_handler_p->hGLRC)) == FALSE)
    { 
        printf("make_current_gl_handler: wglMakeCurrent failed\n");
        ok = 0;
    }
    if(ReleaseDC(gl_handler_p->hWnd, hDC) == 0) 
        printf("make_current_gl_handler: ReleaseDC failed\n");
#else
    glXWaitGL();
    glXWaitX();
    //glGetError(); /* remove existing errors */
    /* wb104 18 Jul 2013 */
    /* on Mac at least if above is not used then last "on" spectrum is not drawn */
    /* that is because glGetError() call below was throwing a 1282 (GL_INVALID_OPERATION) error */
    /* but above glGetError() was crashing some Linux at least */
    ok = glXMakeCurrent(gl_handler_p->display, Tk_WindowId(gl_handler_p->tk_win),
						gl_handler_p->context);
#endif
/*
    glDrawBuffer(GL_BACK);
*/
    //if (!ok || !(glGetError() == GL_NO_ERROR))
    //  return CCPN_FALSE;
    if (!ok)
    {
      printf("not ok\n");
      return CCPN_FALSE;
    }
/*  wb104 17 Mar 2014 */
/*  below often seems to be triggered on Mac (see above) but leaving it out seems to cause no harm */
/*
    GLenum err = glGetError();
    if (err != GL_NO_ERROR)
      return CCPN_FALSE;
*/

    return CCPN_TRUE;
}

void swap_buffers_gl_handler(Gl_handler gl_handler)
{
#ifdef WIN64
    HDC hDC;
#endif
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

/*
    printf("swap_buffers_gl_handler\n");
*/
    if (!make_current_gl_handler(gl_handler))
        return;

#ifdef WIN64
    if((hDC = GetDC(gl_handler_p->hWnd)) == NULL) 
        printf("swap_buffers_gl_handler: GetDC failed\n");
    /* SwapBuffers(hDC); */
    if( wglSwapLayerBuffers(hDC, WGL_SWAP_MAIN_PLANE) == FALSE)
        printf("swap_buffers_gl_handler: wglSwapLayerBuffers failed\n");
    if(ReleaseDC(gl_handler_p->hWnd, hDC) == 0)
        printf("swap_buffers_gl_handler: ReleaseDC failed\n");
#else
    glXSwapBuffers(gl_handler_p->display, Tk_WindowId(gl_handler_p->tk_win));
#endif
/*
    XDrawString(gl_handler->display, Tk_WindowId(gl_handler->tk_win),
			gl_handler->gc, 50, 50, "hello", 5);
*/
}

void test_gl_handler(Gl_handler gl_handler)
{
    int i, n = 100;
    float s = 0.5, dx = 0.5, dy = 0.5, x, y, ang;
    float r = 1.0, g = 0.0, b = 0.0;
/*
GLdouble eqn[4];
*/

/*
    printf("test_gl_handler\n");
*/
/*
    if (!make_current_gl_handler(gl_handler))
      return;
*/

/*
    glEnable(GL_CLIP_PLANE0);
    glGetClipPlane(GL_CLIP_PLANE0, eqn);
    printf("test clip eqn = %3.2f, %3.2f, %3.2f, %3.2f\n",
				eqn[0], eqn[1], eqn[2], eqn[3]);
    printf("clipping is enabled = %d\n", glIsEnabled(GL_CLIP_PLANE0));
*/
    glBegin(GL_LINES);

    for (i = 0; i < n; i++)
    {
	g = (float) i / n;
	glColor3f(r, g, b);
	glVertex2f(dx, dy);
	ang = 2.0 * i * 3.14159265 / n;
	x = s * cos((double) ang) + dx;
	y = s * sin((double) ang) + dy;
	glVertex2f(x, y);
    }

    glEnd();

/*
    XDrawString(gl_handler->display, Tk_WindowId(gl_handler->tk_win),
			gl_handler->gc, 50, 50, "hello", 5);
*/
}

void draw_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1)
{
/*
    glColor3f(0.0, 0.0, 1.0);
    printf("draw_box_gl_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/

/*
    glBegin(GL_LINE_LOOP);
    glVertex2f(x0, y0);
    glVertex2f(x1, y0);
    glVertex2f(x1, y1);
    glVertex2f(x0, y1);
    glEnd();
*/
    draw_clipped_line_gl_handler(gl_handler, x0, y0, x1, y0);
    draw_clipped_line_gl_handler(gl_handler, x1, y0, x1, y1);
    draw_clipped_line_gl_handler(gl_handler, x1, y1, x0, y1);
    draw_clipped_line_gl_handler(gl_handler, x0, y1, x0, y0);
}

void draw_dash_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1)
{
    GLint factor = 1;
    GLushort pattern = determinePattern(2, 6);

/*
    printf("draw_dash_box_gl_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(factor, pattern);

/*
    glBegin(GL_LINE_LOOP);
    glVertex2f(x0, y0);
    glVertex2f(x1, y0);
    glVertex2f(x1, y1);
    glVertex2f(x0, y1);
    glEnd();
*/
    draw_clipped_line_gl_handler(gl_handler, x0, y0, x1, y0);
    draw_clipped_line_gl_handler(gl_handler, x1, y0, x1, y1);
    draw_clipped_line_gl_handler(gl_handler, x1, y1, x0, y1);
    draw_clipped_line_gl_handler(gl_handler, x0, y1, x0, y0);

    glDisable(GL_LINE_STIPPLE);
}

void draw_xor_box_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1)
{
/*
    printf("draw_xor_box_gl_handler: %4.3f, %4.3f, %4.3f, %4.3f\n", x0, y0, x1, y1);
*/
    // NOTE:ED changed to alternating line stipple
    static GLubyte stipple[] = {
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
        0x88, 0x88, 0x88, 0x88, 0x00, 0x00, 0x00, 0x00, 0x22, 0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00,
    };

/*
    start_xor_gl_handler(gl_handler);
*/
    glEnable(GL_POLYGON_STIPPLE);
    glPolygonStipple(stipple);
/*
    draw_box_gl_handler(gl_handler, x0, y0, x1, y1);
*/
    glRectf(x0, y0, x1, y1);
    glDisable(GL_POLYGON_STIPPLE);

    draw_clipped_line_gl_handler(gl_handler, x0, y0, x1, y0);
    draw_clipped_line_gl_handler(gl_handler, x1, y0, x1, y1);
    draw_clipped_line_gl_handler(gl_handler, x1, y1, x0, y1);
    draw_clipped_line_gl_handler(gl_handler, x0, y1, x0, y0);

/*
    finish_xor_gl_handler(gl_handler);
*/
}

void start_xor_gl_handler(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLuint main_call_list = gl_handler_p->main_call_list;
    GLuint xor_call_list = gl_handler_p->xor_call_list;

/*
    printf("start_xor_gl_handler\n");
*/
/*
    clear_xor_gl_handler(gl_handler);
*/

    glClear(GL_COLOR_BUFFER_BIT);
    glCallList(main_call_list);

    glNewList(xor_call_list, GL_COMPILE_AND_EXECUTE);
/* not clear whether glLogicOp has to be before or after glEnable */
/*
    glLogicOp(GL_XOR);
*/
    glEnable(GL_COLOR_LOGIC_OP);
    glLogicOp(GL_XOR);
/*
    glDrawBuffer(GL_FRONT);
*/
}

void finish_xor_gl_handler(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLuint call_list = gl_handler_p->xor_call_list;

/*
    printf("finish_xor_gl_handler\n");
*/
/*
    glDrawBuffer(GL_BACK);
*/

    glEndList();
    glDisable(GL_COLOR_LOGIC_OP);
    if (gl_handler_p->is_double_buffer)
	swap_buffers_gl_handler(gl_handler);
}

void reset_xor_gl_handler(Gl_handler gl_handler)
{
/*
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLuint call_list = gl_handler_p->xor_call_list;
*/

/*
    printf("reset_xor_gl_handler\n");
*/
/*
    glCallList(call_list);
*/
}

void clear_xor_gl_handler(Gl_handler gl_handler)
{
/*
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLuint call_list = gl_handler_p->xor_call_list;
    GLsizei width = gl_handler_p->width;
    GLsizei height = gl_handler_p->height;
*/

/*
    printf("clear_xor_gl_handler\n");
*/

/*
    glCallList(call_list);
    glNewList(call_list, GL_COMPILE);
    glEndList();
*/
}

void start_front_gl_handler(Gl_handler gl_handler)
{
/*
    printf("start_front_gl_handler\n");
*/

    glDrawBuffer(GL_FRONT);
}

void start_back_gl_handler(Gl_handler gl_handler)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GLuint call_list = gl_handler_p->main_call_list;

/*
    printf("start_back_gl_handler\n");
*/

    glDrawBuffer(GL_BACK);

    glNewList(call_list, GL_COMPILE_AND_EXECUTE);
}

void end_back_gl_handler(Gl_handler gl_handler)
{
/*
    printf("end_back_gl_handler\n");
*/

    glEndList();
}

void set_background_gl_handler(Gl_handler gl_handler, float *background)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    COPY_VECTOR(gl_handler_p->background, background, NCOLORS);
}

void draw_text_gl_handler(Gl_handler gl_handler, CcpnString text, float x, float y,
							float a, float b)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GlutFont *glutFont = gl_handler_p->glutFont;
    char *c;
    GLfloat xx, yy;
    GLdouble s, sx, sy;
    static GLfloat fontScale = 119.05 + 33.33;
    float w, h;

/*
    printf("draw_text_gl_handler: %s\n", text);
*/
    if (!glutFont)
	return;

    get_text_size_gl_handler(gl_handler, text, &w, &h);

/* TBD: is this correct?? */
    xx = x - a * w / (gl_handler_p->sx * gl_handler_p->width);
    yy = y - b * h / (gl_handler_p->sy * gl_handler_p->height);

/* not sure if need to push/pop */
    glPushMatrix();
    if (glutFont->isStrokeFont)
    {
        glTranslatef(xx, yy, 0.0);
        s = glutFont->size / fontScale;
	sx = s / (gl_handler_p->sx * gl_handler_p->width);
	sy = s / (gl_handler_p->sy * gl_handler_p->height);
        glScaled(sx, sy, 1.0);
        for (c = text; *c; c++)
	    glutStrokeCharacter(glutFont->font, *c);
    }
    else
    {
        glRasterPos2f(xx, yy);
        for (c = text; *c; c++)
	    glutBitmapCharacter(glutFont->font, *c);
    }
    glPopMatrix();
}

void draw_line_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1)
{
    GLfloat a0 = x0, a1 = x1, b0 = y0, b1 = y1;

    glBegin(GL_LINES);
    glVertex2f(a0, b0);
    glVertex2f(a1, b1);
    glEnd();
}

void draw_clipped_line_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    float xmin = gl_handler_p->x0;
    float ymin = gl_handler_p->y0;
    float xmax = gl_handler_p->x1;
    float ymax = gl_handler_p->y1;

    draw_clipped_line(x0, y0, x1, y1, &drawing_funcs,
		(Generic_ptr) gl_handler, xmin, ymin, xmax, ymax);
}

#define  NCIRCLE_PTS  20
#define  DEGREES_PER_PT  (360 / NCIRCLE_PTS)
#define  RADIANS_PER_PT  (DEGREES_PER_PT * RADIAN_SCALE)

void fill_circle_gl_handler(Gl_handler gl_handler, float x, float y, float r)
{
    fill_ellipse_gl_handler(gl_handler, x, y, r, r);
}

void fill_ellipse_gl_handler(Gl_handler gl_handler, float x, float y, float rx, float ry)
{
    int i;
    float x1, y1, x2, y2, angle;

/*
    printf("fill_ellipse_gl_handler: %2.1f %2.1f %2.1f %2.1f\n", x, y, rx, ry);
*/
    x1 = x + rx;
    y1 = y;
    glBegin(GL_TRIANGLES);
    for (i = 1; i <= NCIRCLE_PTS; i++)
    {
        angle = i * RADIANS_PER_PT;
        x2 = x + rx * (float) cos((double) angle);
        y2 = y + ry * (float) sin((double) angle);
        glVertex2f(x, y);
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
        x1 = x2;
        y1 = y2;
    }

    glEnd();
}

void draw_circle_gl_handler(Gl_handler gl_handler, float x, float y, float r)
{
    draw_ellipse_gl_handler(gl_handler, x, y, r, r);
}

void draw_ellipse_gl_handler(Gl_handler gl_handler, float x, float y, float rx, float ry)
{
    int i;
    float x1, y1, x2, y2, angle;

/*
    printf("draw_ellipse_gl_handler: %2.1f %2.1f %2.1f\n", x, y, r);
*/
    x1 = x + rx;
    y1 = y;
    glBegin(GL_LINE_LOOP);
    for (i = 1; i < NCIRCLE_PTS; i++)
    {
        angle = i * RADIANS_PER_PT;
        x2 = x + rx * (float) cos((double) angle);
        y2 = y + ry * (float) sin((double) angle);
        glVertex2f(x1, y1);
        glVertex2f(x2, y2);
        x1 = x2;
        y1 = y2;
    }

    glEnd();
}

void draw_polyline_gl_handler(Gl_handler gl_handler, Poly_line polyline)
{
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
    GLfloat a, b;

/*
printf("gl_draw_polyline: n = %2d\n", n);
*/

    if (polyline->closed)
	glBegin(GL_LINE_LOOP);
    else
	glBegin(GL_LINE_STRIP);

fflush(stdout);
    for (i = 0; i < n; i++)
    {
	a = v[i].x;
	b = v[i].y;
	glVertex2f(a, b);
    }

    glEnd();
}

void draw_clipped_polyline_gl_handler(Gl_handler gl_handler, Poly_line polyline)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    float xmin = gl_handler_p->x0;
    float ymin = gl_handler_p->y0;
    float xmax = gl_handler_p->x1;
    float ymax = gl_handler_p->y1;

    draw_clipped_polyline(polyline, &drawing_funcs,
		(Generic_ptr) gl_handler, xmin, ymin, xmax, ymax);
}

void draw_dash_line_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
    GLint factor = 1;
    GLushort pattern = determinePattern(dash_length, gap_length);

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(factor, pattern);

    draw_line_gl_handler(gl_handler, x0, y0, x1, y1);

    glDisable(GL_LINE_STIPPLE);
}

void draw_clipped_dash_line_gl_handler(Gl_handler gl_handler,
			float x0, float y0, float x1, float y1,
			int dash_length, int gap_length)
{
    GLint factor = 1;
    GLushort pattern = determinePattern(dash_length, gap_length);

    glEnable(GL_LINE_STIPPLE);
    glLineStipple(factor, pattern);

    draw_clipped_line_gl_handler(gl_handler, x0, y0, x1, y1);

    glDisable(GL_LINE_STIPPLE);
}

void set_color_gl_handler(Gl_handler gl_handler, float *color)
{
    GLfloat red = color[0], green = color[1], blue = color[2];

    glColor3f(red, green, blue);
}

void set_black_gl_handler(Gl_handler gl_handler)
{
    float color[NCOLORS];

    ZERO_VECTOR(color, NCOLORS);

    set_color_gl_handler(gl_handler, color);
}

void set_white_gl_handler(Gl_handler gl_handler)
{
    int i;
    float color[NCOLORS];

    for (i = 0; i < NCOLORS; i++)
        color[i] = 1;

    set_color_gl_handler(gl_handler, color);
}

void set_line_width_gl_handler(Gl_handler gl_handler, float line_width)
{
    GLfloat width = line_width;

    glLineWidth(width);
}

void reset_line_width_gl_handler(Gl_handler gl_handler)
{
    set_line_width_gl_handler(gl_handler, DEFAULT_LINEWIDTH);
}

CcpnStatus has_font_gl_handler(Gl_handler gl_handler, CcpnString name, int *size)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    if (gl_handler_p->glutFont != NULL) {
        *size = gl_handler_p->glutFont->size;
        sprintf(name, gl_handler_p->glutFont->name);
	    return CCPN_TRUE;
    }
	else
	    return CCPN_FALSE;
}

CcpnStatus set_font_gl_handler(Gl_handler gl_handler, CcpnString name, int size)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GlutFont *glutFont;
    int i, imin, d, dmin;

/* look for exact match */
    for (i = 0; i < nbitmapFonts; i++)
    {
	glutFont = glutBitmapFonts + i;
	if ((size == glutFont->size) && equal_strings(name, glutFont->name))
	{
	    gl_handler_p->glutFont = glutFont;
	    return CCPN_OK;
	}
    }

    for (i = 0; i < nstrokeFonts; i++)
    {
	glutFont = glutStrokeFonts + i;
	if (equal_strings(name, glutFont->name))
	{
	    gl_handler_p->glutFont = glutFont;
	    glutFont->size = size;
	    return CCPN_OK;
	}
    }

/* look for size match */
    for (i = 0; i < nbitmapFonts; i++)
    {
	glutFont = glutBitmapFonts + i;
	if (size == glutFont->size)
	{
	    gl_handler_p->glutFont = glutFont;
	    return CCPN_OK;
	}
    }

    if (nstrokeFonts > 0)
    {
	glutFont = glutStrokeFonts;
	gl_handler_p->glutFont = glutFont;
	return CCPN_OK;
    }

/* look for closest size match in font */
    imin = -1;
    for (i = 0; i < nbitmapFonts; i++)
    {
	glutFont = glutBitmapFonts + i;
	if (equal_strings(name, glutFont->name))
	{
	    d = ABS(size - glutFont->size);
	    if ((imin < 0) || (d < dmin))
	    {
		imin = i;
		dmin = d;
	    }
	}
    }

    if (imin >= 0)
    {
	glutFont = glutBitmapFonts + imin;
	gl_handler_p->glutFont = glutFont;
	return CCPN_OK;
    }

/* look for closest size match */
    for (i = 0; i < nbitmapFonts; i++)
    {
	glutFont = glutBitmapFonts + i;
	d = ABS(size - glutFont->size);
	if ((imin < 0) || (d < dmin))
	{
	    imin = i;
	    dmin = d;
	}
    }

    if (imin >= 0) /* will be true if at least one font */
    {
	glutFont = glutBitmapFonts + imin;
	gl_handler_p->glutFont = glutFont;
	return CCPN_OK;
    }

    return CCPN_ERROR;
}

void get_region_gl_handler(Gl_handler gl_handler, float *x0, float *y0,
                                                    float *x1, float *y1)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    *x0 = gl_handler_p->x0;
    *y0 = gl_handler_p->y0;
    *x1 = gl_handler_p->x1;
    *y1 = gl_handler_p->y1;
}

void get_text_size_gl_handler(Gl_handler gl_handler, CcpnString text,
						float *p_w, float *p_h)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;
    GlutFont *glutFont = gl_handler_p->glutFont;
    char *c;
    int w;

    if (!glutFont)
    {
	/* arbitrary */
	*p_h = 10;
	*p_w = 7 * strlen(text);
	return;
    }

    w = 0;
    for (c = text; *c; c++)
	w += glutBitmapWidth(glutFont->font, *c);

    *p_w = w;
    *p_h = glutFont->size;
}

void fill_triangle_gl_handler(Gl_handler gl_handler, float x0, float y0,
				float x1, float y1, float x2, float y2)
{
    glBegin(GL_TRIANGLES);
    glVertex2f(x0, y0);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}


void set_is_double_buffer_gl_handler(Gl_handler gl_handler, Bool is_double_buffer)
{
    Gl_handler_p gl_handler_p = (Gl_handler_p) gl_handler;

    gl_handler_p->is_double_buffer = is_double_buffer;
}

#endif /* IGNORE_GL */


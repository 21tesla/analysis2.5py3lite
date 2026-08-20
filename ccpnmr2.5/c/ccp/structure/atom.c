
/*
======================COPYRIGHT/LICENSE START==========================

atom.c: Part of the CcpNmr Analysis program

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
#include "atom.h"

#include "bond.h"
#include "color.h"

#define  ALLOC_INCR  4

Atom new_atom(float size, char *symbol, char *annotation,
                                        float *x, float *color)
{
    Atom atom;

    MALLOC_NEW(atom, struct Atom, 1);

    atom->size = size;
    atom->symbol = NULL;
    atom->annotation = NULL;
    atom->bonds = NULL;

    COPY_VECTOR(atom->x, x, ATOM_NDIMS);
    set_color_atom(atom, color);

    if (set_symbol_atom(atom, symbol) == CCPN_ERROR)
       return NULL;

    if (set_annotation_atom(atom, annotation) == CCPN_ERROR)
       return NULL;

    atom->nbonds = 0;
    MALLOC_NEW(atom->bonds, Bond, ALLOC_INCR);
    atom->nbonds_alloc = ALLOC_INCR;

    atom->is_drawn = CCPN_TRUE;
    atom->have_annotation_color = CCPN_FALSE;

    return atom;
}
 
void delete_atom(Atom atom)
{
/*
    printf("delete_atom\n");
*/

    if (atom)
    {
        FREE(atom->bonds, Bond);
        FREE(atom->symbol, char);
        FREE(atom->annotation, char);
        FREE(atom, struct Atom);
    }
}

void set_size_atom(Atom atom, float size)
{
    atom->size = size;
}

CcpnStatus set_symbol_atom(Atom atom, char *symbol)
{
    FREE(atom->symbol, char);

    STRING_MALLOC_COPY(atom->symbol, symbol);

    return CCPN_OK;
}

CcpnStatus set_annotation_atom(Atom atom, char *annotation)
{
    FREE(atom->annotation, char);

    STRING_MALLOC_COPY(atom->annotation, annotation);

    return CCPN_OK;
}

void set_color_atom(Atom atom, float *color)
{
    COPY_VECTOR(atom->color, color, ATOM_NCOLORS);
}

void turn_on_atom(Atom atom)
{
    atom->is_drawn = CCPN_TRUE;
}

void turn_off_atom(Atom atom)
{
    atom->is_drawn = CCPN_FALSE;
}

void set_annotation_color_atom(Atom atom, float *color)
{
    if (color)
    {
        atom->have_annotation_color = CCPN_TRUE;
        COPY_VECTOR(atom->annotation_color, color, ATOM_NCOLORS);
    }
    else
    {
      atom->have_annotation_color = CCPN_FALSE;
    }
}

CcpnStatus add_bond_atom(Atom atom, Bond bond)
{
    int nalloc = atom->nbonds_alloc, n = atom->nbonds;

    if (n == nalloc)
    {
        nalloc += ALLOC_INCR;
        REALLOC(atom->bonds, Bond, nalloc);
        atom->nbonds_alloc = nalloc;
    }

    atom->bonds[n] = bond;
    atom->nbonds++;

    return CCPN_OK;
}

CcpnStatus remove_bond_atom(Atom atom, Bond bond)
{
    int i, n = atom->nbonds - 1;

    /* assume that later bonds more likely to be removed, so work from back */
    for (i = n; i >= 0; i--)
    {
        if (atom->bonds[i] == bond)
            break;
    }

    if (i < 0)
        return CCPN_ERROR;

    /* just set last bond to fill this slot since order does not matter */
    atom->bonds[i] = atom->bonds[n];
    atom->nbonds--;

    return CCPN_OK;
}

#define  MAX_Z  10.0
#define  MIN_PARAM  0.2
#define  MAX_PARAM  1.0
#define  A  (0.5*(MAX_PARAM+MIN_PARAM))
#define  B  (0.5*(MAX_PARAM-MIN_PARAM)/MAX_Z)

/* param = MAX_PARAM for z >= MAX_Z
 * param = MIN_PARAM for z <= - MAX_Z
 * and linearly interpolated in between
 * color = param*original_color + (1-param)*background_color
 * so param = 1 means ignore background
 * and param = 0 means just use background
 */
static float get_depth_param_atom(Atom atom, float depth)
{
    float z = atom->x[2] - depth;

    z = MAX(-MAX_Z, MIN(MAX_Z, z));

    return (A + B*z);
}

void draw_atom(Atom atom, float camera, float depth,
                Drawing_funcs *drawing_funcs, Generic_ptr data)
{
    int i;
    float x = atom->x[0], y = atom->x[1], z = atom->x[2], r = atom->size;
    float j, param, param2, r2, smooth = 0.005;
    float fill_color[ATOM_NCOLORS];
    float fill_color2[ATOM_NCOLORS];
    float arc_color[ATOM_NCOLORS];
    float background_color[ATOM_NCOLORS], symbol_color[ATOM_NCOLORS];
    float field_depth = -4.0;

    if (!atom->is_drawn)
       return;

    /* Perspective transform */
    z = z-camera;
    
    if (z >= 0) {
      return;
    }
    
    x = field_depth *x/z;
    y = field_depth *y/z;
    r = field_depth *r/z;
    /* End perspective transform */

    if (r > smooth)
    {
        (drawing_funcs->get_background)(data, background_color);

        param = get_depth_param_atom(atom, depth);

        for (i = 0; i < ATOM_NCOLORS; i++)
            fill_color[i] = param * atom->color[i] + (1.0-param) * background_color[i];
         
        for (j=0.0; j<r; j+= smooth) {
          param = j/r;
          r2 = param*0.20*r;
          
          if (param < 0.5) {
            param2 = 0.25 + 1.5*param;
            for (i = 0; i < ATOM_NCOLORS; i++)
                fill_color2[i] = param2 * fill_color[i] + (1.0-param2) * background_color[i];
         
          }
          else {
            param2 = 0.7*(param-0.5);
            for (i = 0; i < ATOM_NCOLORS; i++)
                fill_color2[i] = (1.0-param2) * fill_color[i] + param2;
          }
           
          (drawing_funcs->set_draw_color)(data, fill_color2);
          (drawing_funcs->fill_circle)(data, x+r2, y+r2, r*cos(param*1.5707963));

        }
        
        /*
        param = 0.75;
        inverted_grey_color(arc_color, fill_color); 
        for (i = 0; i < ATOM_NCOLORS; i++)
            arc_color[i] = param * fill_color[i] + (1.0-param) *  background_color[i];
        (drawing_funcs->set_draw_color)(data, arc_color);
        (drawing_funcs->draw_circle)(data, x, y, r);
      */
    }

    if (atom->annotation && *(atom->annotation))
    {
        if (atom->have_annotation_color)
            (drawing_funcs->set_draw_color)(data, atom->annotation_color);
        (drawing_funcs->draw_text)(data, atom->annotation, x+r, y-r, 0.0, 0.0);
    }

    if (atom->symbol && *(atom->symbol))
    {
        inverted_grey_color(symbol_color, fill_color);
/*
        for (i = 0; i < ATOM_NCOLORS; i++)
            symbol_color[i] = param * symbol_color[i] + (1.0-param) * background_color[i];
*/
        (drawing_funcs->set_draw_color)(data, symbol_color);
        (drawing_funcs->draw_text)(data, atom->symbol, x, y, 0.4, 0.5);
    }
}

void translate_atom(Atom atom, float *delta)
{
    ADD_VECTORS(atom->x, atom->x, delta, ATOM_NDIMS);
}

void rotate_atom(Atom atom, float **matrix, float *origin)
{
    int i, j;
    float x[ATOM_NDIMS];

    for (i = 0; i < ATOM_NDIMS; i++)
    {
        x[i] = origin[i];
        for (j = 0; j < ATOM_NDIMS; j++)
           x[i] += matrix[i][j] * (atom->x[j] - origin[j]);
    }

    COPY_VECTOR(atom->x, x, ATOM_NDIMS);
}

void zoom_atom(Atom atom, float scale)
{

    if (scale == 0)
        return;
        
    SCALE_VECTOR(atom->x, atom->x, scale, ATOM_NDIMS);
    
}

void set_coords_atom(Atom atom, float *coords)
{
    int i;

    for (i = 0; i < ATOM_NDIMS; i++)
	atom->x[i] = coords[i];
}

Bool within_xy_tol_atom(Atom atom, float x, float y, float tol)
{
    float dx = atom->x[0] - x;
    float dy = atom->x[1] - y;

/*
    printf("dx = %f, dy = %f, tol = %f\n", dx, dy, tol);
*/

    return (dx*dx+dy*dy < tol*tol) ? CCPN_TRUE : CCPN_FALSE;
}


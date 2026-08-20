
/*
======================COPYRIGHT/LICENSE START==========================

bond.c: Part of the CcpNmr Analysis program

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
#include "bond.h"

Bond new_bond(Atom atom1, Atom atom2, float *color)
{
    Bond bond;

    MALLOC_NEW(bond, struct Bond, 1);

    bond->atom1 = atom1;
    bond->atom2 = atom2;

    set_color_bond(bond, color);
    bond->line_width = -1.0; /* anything negative */
    bond->line_style = NORMAL_LINE_STYLE;
    bond->annotation = NULL;

    return bond;
}
 
void delete_bond(Bond bond)
{
/*
    printf("delete_bond\n");
*/

    if (bond)
    {
        FREE(bond->annotation, char);
        FREE(bond, struct Bond);
    }
}

void set_color_bond(Bond bond, float *color)
{
    if (color)
    {
        bond->have_color = CCPN_TRUE;
        COPY_VECTOR(bond->color, color, BOND_NCOLORS);
    }
    else
    {
        bond->have_color = CCPN_FALSE;
    }
}

void set_line_width_bond(Bond bond, float line_width)
{
    bond->line_width = line_width;
}

void set_line_style_bond(Bond bond, int line_style)
{
    if ((line_style >= 0) && (line_style < NLINE_STYLES))
        bond->line_style = line_style;
}

CcpnStatus set_annotation_bond(Bond bond, char *annotation)
{
    FREE(bond->annotation, char);

    STRING_MALLOC_COPY(bond->annotation, annotation);

    return CCPN_OK;
}

Atom get_other_atom_bond(Bond bond, Atom atom)
{
    if (bond->atom1 == atom)
        return bond->atom2;
    else if (bond->atom2 == atom)
        return bond->atom1;
    else
        return NULL;
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
static float get_depth_param_bond(Bond bond, float depth)
{
    float z = 0.5 * (bond->atom1->x[2] + bond->atom2->x[2]) - depth;

    z = MAX(-MAX_Z, MIN(MAX_Z, z));

    return (A + B*z);
}

void draw_bond(Bond bond, float camera, float depth, Drawing_funcs *drawing_funcs,
                                                        Generic_ptr data)
{
    Atom atom1 = bond->atom1, atom2 = bond->atom2;
    float a1 = atom1->x[0], b1 = atom1->x[1]; 
    float a2 = atom2->x[0], b2 = atom2->x[1]; 
    float c1 = atom1->x[2], c2 = atom2->x[2];
    float a3, b3, a = 0.0, b = 0.0, r, param, line_width;
    float color[BOND_NCOLORS], background_color[BOND_NCOLORS];
    int i;
    float field_depth = -4.0;

    /* Perspective transform */
    c1 = c1-camera;
    c2 = c2-camera;
    
    if (c1 >= 0) {
      return;
    }
    if (c2 >=0) {
      return;
    }
    
    a1 = field_depth * a1/c1; 
    b1 = field_depth * b1/c1;

    a2 = field_depth * a2/c2;
    b2 = field_depth * b2/c2;
    
    line_width = 8.0 * (float) bond->line_width;
    line_width = 2.0 * field_depth * line_width/(c1+c2);
    
    /* End perspective transform */

    if (bond->line_width >= 0)
        (drawing_funcs->set_line_width)(data, 1+(int) line_width);

    if (bond->line_style != NORMAL_LINE_STYLE)
        (drawing_funcs->set_line_style)(data, bond->line_style);

    (drawing_funcs->get_background)(data, background_color);

    param = get_depth_param_bond(bond, depth);
    if (bond->have_color)
    {
        for (i = 0; i < BOND_NCOLORS; i++)
            color[i] = param * bond->color[i] + (1.0-param) * background_color[i];
        (drawing_funcs->set_draw_color)(data, color);
        (drawing_funcs->draw_line)(data, a1, b1, a2, b2);
    }
    else
    {
        a3 = 0.5 * (a1 + a2);
        b3 = 0.5 * (b1 + b2);

        for (i = 0; i < BOND_NCOLORS; i++)
            color[i] = param * atom1->color[i] + (1.0-param) * background_color[i];
        (drawing_funcs->set_draw_color)(data, color);
        (drawing_funcs->draw_line)(data, a1, b1, a3, b3);

        for (i = 0; i < BOND_NCOLORS; i++)
            color[i] = param * atom2->color[i] + (1.0-param) * background_color[i];
        (drawing_funcs->set_draw_color)(data, color);
        (drawing_funcs->draw_line)(data, a3, b3, a2, b2);
    }

    if (bond->annotation && *(bond->annotation))
    {
        if (bond->have_color)
        {
            a3 = 0.5 * (a1 + a2);
            b3 = 0.5 * (b1 + b2);
        }

        /* 
        r = 0.6 * 0.5 * (bond->atom1->size + bond->atom2->size)
      */

        /* Perspective transform */
        r = field_depth * 0.6 * (atom1->size + atom2->size) / (c1+c2);
        a3 += r;

        if (a1 < a2)
        {
            if (b1 < b2)
                r = -r;
        }
        else
        {
            if (b1 > b2)
                r = -r;
        }

        b3 += r;
        (drawing_funcs->draw_text)(data, bond->annotation, a3, b3, a, b);
    }

    if (bond->line_style != NORMAL_LINE_STYLE)
        (drawing_funcs->set_line_style)(data, NORMAL_LINE_STYLE);

    if (bond->line_width >= 0)
        (drawing_funcs->set_line_width)(data, DEFAULT_BOND_WIDTH);
}

Bool within_xy_tol_bond(Bond bond, float x, float y, float tol, float camera, float *z)
{
    Atom atom1 = bond->atom1;
    Atom atom2 = bond->atom2;
    float x1, y1, z1, x2, y2, z2, z1t, z2t, dx, dy;
    float ax, ay, px, py, b, c, vx, vy, det, lambda;
    float field_depth = -4.0;

#define SMALL_D  1.0e-3


    z1 = atom1->x[2];
    z2 = atom2->x[2];

    /* Perspective transforrm */
    if (z1 >= camera) {
      return CCPN_FALSE;
    }  
    if (z2 >= camera) {
      return CCPN_FALSE;    
    }
    z1t = z1-camera;
    z2t = z2-camera;
    
    x1 = field_depth * atom1->x[0]/z1t;
    y1 = field_depth * atom1->x[1]/z1t;
    x2 = field_depth * atom2->x[0]/z2t;
    y2 = field_depth * atom2->x[1]/z2t;
    
    /* End perspective transform */

    dx = x2 - x1;
    dy = y2 - y1;

    if ((ABS(dx) < SMALL_D) && (ABS(dy) < SMALL_D))
    {
	/* points are on top of each other */
	vx = x1;
	vy = y1;
	*z = MAX(z1, z2);
    }
    else
    {
	if (ABS(dx) > ABS(dy))
	{
	    ax = -dy/dx;
	    py = -ax;
	    px = ay = 1;
	}
	else
	{
	    py = ax = 1;
	    ay = -dx/dy;
	    px = -ay;
	}

	b = ax*x1 + ay*y1;
	c = ax*x + ay*y;

	det = ax*py - ay*px;
	vx = (py*b - ay*c) / det;
	vy = (-px*b + ax*c) / det;

	lambda = (px*(vx-x1) + py*(vy-y1)) / (px*dx + py*dy);

	if (lambda < 0)
	{
	    vx = x1;
	    vy = y1;
	    *z = z1;
	}
	else if (lambda > 1)
	{
	    vx = x2;
	    vy = y2;
	    *z = z2;
	}
	else
	{
	    *z = (1-lambda)*z1 + lambda*z2;
	}
    }

    dx = x - vx;
    dy = y - vy;

/*
    printf("dx = %f, dy = %f, tol = %f\n", dx, dy, tol);
*/

    return (dx*dx+dy*dy < tol*tol) ? CCPN_TRUE : CCPN_FALSE;
}


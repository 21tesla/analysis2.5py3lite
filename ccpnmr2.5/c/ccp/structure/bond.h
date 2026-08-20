/*
======================COPYRIGHT/LICENSE START==========================

bond.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_bond
#define _incl_bond

#include "macros.h"
#include "types.h"

#include "atom.h"
#include "drawing_funcs.h"

#define BOND_NCOLORS  3

#define DEFAULT_BOND_WIDTH  3.0

typedef struct Bond
{
    Atom atom1;
    Atom atom2;
    Bool have_color;
    float color[BOND_NCOLORS];
    float line_width;
    int line_style;
    char *annotation;
}   *Bond;

extern Bond new_bond(Atom atom1, Atom atom2, float *color);
                    
extern void delete_bond(Bond bond);

extern void set_color_bond(Bond bond, float *color);

extern void set_line_width_bond(Bond bond, float line_width);

extern void set_line_style_bond(Bond bond, int line_style);

extern CcpnStatus set_annotation_bond(Bond bond, char *annotation);

extern Atom get_other_atom_bond(Bond bond, Atom atom);

extern void draw_bond(Bond bond, float camera, float depth,
                         Drawing_funcs *drawing_funcs, Generic_ptr data);

extern Bool within_xy_tol_bond(Bond bond, float x, float y, float tol,
                                 float camera,  float *z);

#endif /* _incl_bond */

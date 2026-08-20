/*
======================COPYRIGHT/LICENSE START==========================

structure.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_structure
#define _incl_structure

#include "macros.h"
#include "types.h"

#include "atom.h"
#include "bond.h"
#include "drawing_funcs.h"

typedef struct Structure
{
    int natoms;
    int natoms_alloc;
    Atom *atoms;
    int nbonds;
    int nbonds_alloc;
    float camera;
    Bond *bonds;
}   *Structure;

extern Structure new_structure(void);

extern void delete_structure(Structure structure);

extern CcpnStatus add_atom_structure(Structure structure, Atom atom);

extern CcpnStatus remove_atom_structure(Structure structure, Atom atom);

extern Atom nearest_atom_structure(Structure structure, float x, float y, float tol);

extern CcpnStatus add_bond_structure(Structure structure, Bond bond);

extern CcpnStatus remove_bond_structure(Structure structure, Bond bond);

extern Bond nearest_bond_structure(Structure structure, float x, float y, float tol);

extern void draw_structure(Structure structure, Drawing_funcs *drawing_funcs, Generic_ptr data);

extern void translate_structure(Structure structure, float *delta);

extern void rotate_structure(Structure structure, float **matrix, float *origin);

extern void zoom_structure(Structure structure, float delta);

extern void move_to_center_structure(Structure structure);

extern void focus_on_atom_structure(Structure structure, Atom atom);

#endif /* _incl_structure */

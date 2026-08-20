/*
======================COPYRIGHT/LICENSE START==========================

atom_coord_list.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_atom_coord_list
#define _incl_atom_coord_list

#include "atom_coord.h"

typedef struct Atom_coord_list
{
    int natom_coords;
    int natom_coords_alloc;
    Atom_coord *atom_coords;
}   *Atom_coord_list;

extern Atom_coord_list new_atom_coord_list(int natoms, float *mass,
					float *x, float *y, float *z);

extern void delete_atom_coord_list(Atom_coord_list atom_coord_list);

extern void clear_atom_coord_list(Atom_coord_list atom_coord_list);

extern CcpnStatus add_atom_coord_list(Atom_coord_list atom_coord_list,
		float mass, float x, float y, float z, CcpnString error_msg);

extern CcpnStatus append_atom_coord_list(Atom_coord_list atom_coord_list,
				Atom_coord atom_coord, CcpnString error_msg);

#endif /* _incl_atom_coord_list */

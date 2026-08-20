
/*
======================COPYRIGHT/LICENSE START==========================

atom_coord_list.c: Part of the CcpNmr Analysis program

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
#include "atom_coord_list.h"

#define  NALLOC  100

Atom_coord_list new_atom_coord_list(int natoms, float *mass,
					float *x, float *y, float *z)
{
    int i;
    Atom_coord atom_coord;
    Atom_coord_list atom_coord_list;

    MALLOC_NEW(atom_coord_list, struct Atom_coord_list, 1);

    if (natoms > 0)
    {
	MALLOC_NEW(atom_coord_list->atom_coords, Atom_coord, natoms);
	atom_coord_list->natom_coords_alloc = natoms;
    }
    else
    {
	atom_coord_list->atom_coords = NULL;
	atom_coord_list->natom_coords_alloc = 0;
    }

    atom_coord_list->natom_coords = natoms;

    for (i = 0; i < atom_coord_list->natom_coords_alloc; i++)
	atom_coord_list->atom_coords[i] = NULL;

    for (i = 0; i < natoms; i++)
    {
	atom_coord = new_atom_coord(mass[i], x[i], y[i], z[i]);
	if (!atom_coord)
	    return NULL;

	atom_coord_list->atom_coords[i] = atom_coord;
    }

    return atom_coord_list;
}
 
void delete_atom_coord_list(Atom_coord_list atom_coord_list)
{
    int i;

    if (atom_coord_list)
    {
	for (i = 0; i < atom_coord_list->natom_coords; i++)
	    delete_atom_coord(atom_coord_list->atom_coords[i]);

	FREE(atom_coord_list->atom_coords, Atom_coord);

	FREE(atom_coord_list, struct Atom_coord_list);
    }
}

void clear_atom_coord_list(Atom_coord_list atom_coord_list)
{
    if (atom_coord_list)
    {
	FREE(atom_coord_list->atom_coords, Atom_coord);

	FREE(atom_coord_list, struct Atom_coord_list);
    }
}

CcpnStatus add_atom_coord_list(Atom_coord_list atom_coord_list,
		float mass, float x, float y, float z, CcpnString error_msg)
{
    int n;
    int natoms = atom_coord_list->natom_coords;
    int nalloc = atom_coord_list->natom_coords_alloc;
    Atom_coord atom_coord;

    if (natoms >= nalloc)
    {
	sprintf(error_msg, "allocating atom coords memory");
	if (natoms == 0)
	{
	    n = NALLOC;
	    MALLOC(atom_coord_list->atom_coords, Atom_coord, n);
	}
	else
	{
	    n = nalloc + NALLOC;
	    REALLOC(atom_coord_list->atom_coords, Atom_coord, n);
	}

	atom_coord_list->natom_coords_alloc = n;
    }

    atom_coord = new_atom_coord(mass, x, y, z);

    if (!atom_coord)
	RETURN_ERROR_MSG("allocating atom coord memory");

    atom_coord_list->atom_coords[natoms] = atom_coord;
    atom_coord_list->natom_coords++;

    return CCPN_OK;
}

CcpnStatus append_atom_coord_list(Atom_coord_list atom_coord_list,
				Atom_coord atom_coord, CcpnString error_msg)
{
    int n;
    int natoms = atom_coord_list->natom_coords;
    int nalloc = atom_coord_list->natom_coords_alloc;

    if (natoms >= nalloc)
    {
	sprintf(error_msg, "allocating atom coords memory");
	if (natoms == 0)
	{
	    n = NALLOC;
	    MALLOC(atom_coord_list->atom_coords, Atom_coord, n);
	}
	else
	{
	    n = nalloc + NALLOC;
	    REALLOC(atom_coord_list->atom_coords, Atom_coord, n);
	}

	atom_coord_list->natom_coords_alloc = n;
    }

    atom_coord_list->atom_coords[natoms] = atom_coord;
    atom_coord_list->natom_coords++;

    return CCPN_OK;
}


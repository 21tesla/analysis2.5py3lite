/*
======================COPYRIGHT/LICENSE START==========================

struct_util.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_struct_util
#define _incl_struct_util

#include "macros.h"
#include "types.h"

#define COORD_NDIMS  3

/* translate coords to center of mass, using weight if defined */
/* input: coords(atom, dim), weight(atom) */
/* output: coords translated */
extern void translate_coordinates(int natoms, float **coords, float *weight);

/* align coords1 to coords2, after translating to center of mass */
/* input: coords1(atom, dim), coords2(atom, dim), weight(atom) */
/* output: rotation(dim, dim), err */
extern CcpnStatus align_translate_coordinates(int natoms, float **coords1, float **coords2,
                         float *weight, float **rotation, float *err, CcpnString error_msg);

/* align coords1 to coords2 and fill in the alignment rotation matrix */
/* no translation is done */
/* if weight is NULL then all coords taken to have same weight */
/* also returns err, the weighted squared error between the alignments */
/* input: coords1(atom, dim), coords2(atom, dim), weight(atom) */
/* output: rotation(dim, dim), err */
extern CcpnStatus align_coordinates(int natoms, float **coords1, float **coords2,
                float *weight, float **rotation, float *err, CcpnString error_msg);

/* align ensemble of coords */
/* unlike align_coordinates, here weight cannot be null */
/* input: coords(ensemble, atom, dim), weight(atom) */
/* output: updated coords(ensemble, atom, dim), atom_rmsd(atom), ensemble_rmsd(ensemble), toterr */
extern CcpnStatus align_ensemble(int nensembles, int natoms, float ***coords,
                        float *weight, float *atom_rmsd, float *ensemble_rmsd,
                        float *toterr, CcpnString error_msg);

#endif /* _incl_struct_util */

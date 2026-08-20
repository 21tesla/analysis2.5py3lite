
/*
======================COPYRIGHT/LICENSE START==========================

contour_data.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_contour_data
#define _incl_contour_data

#include "block_file.h"
#include "contour_levels.h"
#include "int_array.h"
#include "polyline.h"
#include "store_file.h"
#include "store_handler.h"

typedef struct Contour_data
{
    Int_array block;
    int nplanes;
    int nlevels;
    int **npolylines; /* indexed by (plane, level) */
    Poly_line ***polylines; /* indexed by (plane, level, poly) */
    int cum_planes[MAX_NDIM]; /* includes plane dimensions */
    int nvertices; /* total, for information only */
    int npolys; /* total, for information only */
}   *Contour_data;

extern Contour_data new_contour_data(int xdim, int ydim,
		Block_file block_file, Store_file store_file, Int_array block,
		Contour_levels contour_levels, int ncomponents, int *components,
		Bool transposed, CcpnString error_msg);

extern void delete_contour_data(Contour_data contour_data);

#endif /* _incl_contour_data */

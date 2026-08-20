
/*
======================COPYRIGHT/LICENSE START==========================

store_file.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_store_file
#define _incl_store_file

#include "macros.h"
#include "types.h"

#include "consts.h"

#include "polyline.h"

typedef struct Store_file
{
    FILE *fp;
    Bool swap;
    int header_size;
    int dir_size;
    int *directory;
    int ndim;
    int xdim;
    int ydim;
    Bool have_pos;
    Bool have_neg;
    int block_size[MAX_NDIM];
    int npoints[MAX_NDIM];
    int first[MAX_NDIM];
    int last[MAX_NDIM];
    int nblocks[MAX_NDIM];
    int block_min[MAX_NDIM];
    int cumul_dir[MAX_NDIM];
}   *Store_file;

/* called function takes over ownership of polylines (even if error) */
typedef CcpnStatus (*Store_poly_func)(void *user_data, int *block, int *plane,
	int level, int npolylines, Poly_line *polylines, CcpnString error_msg);

extern Store_file new_store_file(CcpnString fileName, int ndim,
	int xdim, int ydim, int *block_size, CcpnString error_msg);

extern void delete_store_file(Store_file store_file);

/* block is file block (not store block, which can be offset because of first) */
/* plane is block plane (not store block plane, which can be offset because of first) */
extern CcpnStatus process_contours_store_file(Store_file store_file,
			int *block, int *plane, Store_poly_func poly_func,
                        void *user_data, Bool transposed, CcpnString error_msg);

#endif /* _incl_store_file */

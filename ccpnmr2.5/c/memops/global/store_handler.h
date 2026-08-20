
/*
======================COPYRIGHT/LICENSE START==========================

store_handler.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_store_handler
#define _incl_store_handler

#include "macros.h"
#include "types.h"

#include "consts.h"
#include "drawing_funcs.h"

typedef struct Store_handler
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
    int cumul_dir[MAX_NDIM];
    /* rest is for working storage */
    int block[MAX_NDIM];
    int plane[MAX_NDIM];
    int first_plane[MAX_NDIM];
    int offset; /* current offset of data on disk relative to end of directory */
    Bool is_pos;  /* current level is >= 0 */
    int npos_polys;
    int nneg_polys;
}   *Store_handler;

extern Drawing_funcs *store_drawing_funcs(void);

extern Store_handler new_store_handler(CcpnString fileName,
					Bool swap, CcpnString error_msg);

extern void delete_store_handler(Store_handler store_handler);

extern void draw_polyline_store_handler(Store_handler store_handler,
					Poly_line polyline);

extern void draw_clipped_polyline_store_handler(Store_handler store_handler,
					Poly_line polyline);

extern CcpnStatus init_store_save(Store_handler store_handler,
			int ndim, int xdim, int ydim,
                        int *npoints, int *first, int *last,
			int *block_size, int *nblocks,
                        int nlevels, float *levels, CcpnString error_msg);

extern void init_store_block(Store_handler store_handler, int *block);

extern void init_store_plane(Store_handler store_handler, int *plane);

extern void init_store_level(Store_handler store_handler, float level);

extern void end_store_plane(Store_handler store_handler);

extern void end_store_save(Store_handler store_handler);

#endif /* _incl_store_handler */

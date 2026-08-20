
/*
======================COPYRIGHT/LICENSE START==========================

contour_file.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_contour_file
#define _incl_contour_file

#include "block_file.h"
#include "contour_data.h"
#include "contour_levels.h"
#include "contour_style.h"
#include "drawing_funcs.h"
#include "hash_table.h"
#include "mem_cache.h"
#include "shape_file.h"
#include "store_file.h"
#include "store_handler.h"

typedef struct Contour_file
{
    int xdim, ydim;
    Block_file block_file;
    Store_file store_file;
    Contour_levels contour_levels;
    int ncomponents; /* if block_file a SHAPE_FILE */
    int *components; /* if block_file a SHAPE_FILE */
    Hash_table contour_table;
    Mem_cache mem_cache;
    Bool transposed; /* for stored files data can be transposed */
    Bool clearing; /* just needed to help when clearing out contour_table */
}   *Contour_file;

extern Contour_file new_contour_file(int xdim, int ydim,
		Block_file block_file, Store_file store_file,
		Mem_cache mem_cache, Bool transposed, CcpnString error_msg);

extern void delete_contour_file(Contour_file contour_file);

/* draws contours in specified region using specified drawing_funcs */
/* assumes first >= 0 and last <= contour_file.block_file.points */
/* draws all contours with first <= position < last */
/* assumes drawing is externally clipped outside this region in (xdim, ydim) */
/* (contours defined by block, and no clipping in (xdim, ydim) in function */
extern CcpnStatus region_contour_file(Contour_file contour_file,
		int *first, int *last, Abort_func abort_func,
		Contour_levels contour_levels, Contour_style contour_style,
		int ncomponents, int *components, Drawing_funcs *drawing_funcs,
		Generic_ptr data, CcpnString error_msg);

#endif /* _incl_contour_file */

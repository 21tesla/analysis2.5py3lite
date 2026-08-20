
/*
======================COPYRIGHT/LICENSE START==========================

slice_file.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_slice_file
#define _incl_slice_file

#include "block_file.h"
#include "drawing_funcs.h"
#include "mem_cache.h"

typedef struct Slice_file
{
    int orient;
    int dim;
    Block_file block_file;
    Mem_cache mem_cache;
}   *Slice_file;

extern Slice_file new_slice_file(int orient, int dim,
		Block_file block_file, Mem_cache mem_cache);

extern void delete_slice_file(Slice_file slice_file);

/* draws slices in specified region using specified drawing_funcs */
/* assumes first >= 0 and last <= slice_file.block_file.points */
/* draws all slices with first <= position < last */
/* assumes drawing is externally clipped outside this region */
extern CcpnStatus draw_slice_file(Slice_file slice_file,
		int first, int last, float *position,
		Drawing_funcs *drawing_funcs, Generic_ptr data,
		CcpnString error_msg);

/* draws all slices in specified region using specified drawing_funcs */
/* assumes first >= 0 and last <= slice_file.block_file.points */
/* draws all slices with first <= position < last */
/* assumes drawing is externally clipped outside this region */
extern CcpnStatus draw_all_slice_file(Slice_file slice_file,
		int *first, int *last,
                int ncomponents, int *components,
		Drawing_funcs *drawing_funcs, Generic_ptr data,
		CcpnString error_msg);

#endif /* _incl_slice_file */

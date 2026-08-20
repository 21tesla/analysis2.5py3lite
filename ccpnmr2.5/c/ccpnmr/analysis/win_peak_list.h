
/*
======================COPYRIGHT/LICENSE START==========================

win_peak_list.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_win_peak_list
#define _incl_win_peak_list

#include "drawing_funcs.h"
#include "peak_list.h"

#define  DRAW_UNIFORM_METHOD	0
#define  DRAW_GLOBAL_METHOD	1
#define  DRAW_PEAKLIST_METHOD	2
#define  DRAW_LINE_WIDTH_METHOD	3
#define  NDRAW_METHODS		4
/*
#define  DRAW_BOX_WIDTH_METHOD	4
#define  NDRAW_METHODS		5
*/

typedef struct Win_peak_list
{
    Peak_list peak_list;
    Bool hasValueAxis;
    Bool isSymbolDrawn;
    Bool isTextDrawn;
    Bool isTextPointerDrawn;
}   *Win_peak_list;

extern Win_peak_list new_win_peak_list(Peak_list peak_list, Bool hasValueAxis);

extern void delete_win_peak_list(Win_peak_list win_peak_list);

extern void is_symbol_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isSymbolDrawn);
 
extern void is_text_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isTextDrawn);

extern void is_text_pointer_drawn_win_peak_list(Win_peak_list win_peak_list, Bool isTextPointerDrawn);

/* draws all peaks with first <= position < last */
extern CcpnStatus draw_win_peak_list(Win_peak_list win_peak_list, int xdim, int ydim,
		float xpix, float ypix, float xscale, float yscale,
		float *first, float *last,
		int draw_method, float intensity_max, float volume_max,
		float *bg_color, float *center, float *thickness, int *tile,
		Drawing_funcs *drawing_funcs, Generic_ptr data,
		CcpnString error_msg);

#endif /* _incl_win_peak_list */

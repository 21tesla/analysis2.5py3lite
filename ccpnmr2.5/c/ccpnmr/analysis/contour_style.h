
/*
======================COPYRIGHT/LICENSE START==========================

contour_style.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_contour_style
#define _incl_contour_style

#include "macros.h"
#include "types.h"

#define  NORMAL_LINE_STYLE	0
#define  DASHED_LINE_STYLE	1
#define  NLINE_STYLES		2

typedef struct Contour_style
{
    int npos_colors;
    float **pos_colors;
    int nneg_colors;
    float **neg_colors;
    int pos_line_style;
    int neg_line_style;
}   *Contour_style;

/* codes takes ownership of memory */
extern Contour_style new_contour_style(int npos_colors, float **pos_colors,
			int nneg_colors, float **neg_colors,
			int pos_line_style, int neg_line_style);

extern void delete_contour_style(Contour_style contour_style);

#endif /* _incl_contour_style */

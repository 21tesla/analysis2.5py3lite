
/*
======================COPYRIGHT/LICENSE START==========================

peak_list.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_peak_list
#define _incl_peak_list

#include "block_file.h"
#include "peak.h"

typedef struct Peak_list
{
    int ndim;
    int *npoints; /* number of points (for aliasing) */
    float *color; /* (r, g, b) in range [0, 1] */
    int symbol;
    Peak *peaks;
    int npeaks;
    int npeaks_alloc;
}   *Peak_list;

typedef struct Diagonal_exclusion
{
    int dim1;
    int dim2;
    float a1;
    float a2;
    float b12;
    float d;
}   *Diagonal_exclusion;

extern Peak_list new_peak_list(int ndim, int *npoints);

extern void delete_peak_list(Peak_list peak_list);

extern void is_selected_peak_list(Peak_list peak_list, Bool isSelected);

/* peak_list takes copy of color */
extern void set_color_peak_list(Peak_list peak_list, float *color);

extern CcpnStatus set_symbol_peak_list(Peak_list peak_list, int symbol,
						CcpnString error_msg);

extern CcpnStatus add_peak_peak_list(Peak_list peak_list,
				Peak *p_peak, CcpnString error_msg);

extern CcpnStatus remove_peak_peak_list(Peak_list peak_list, Peak peak,
						CcpnString error_msg);

extern void remove_selected_peak_list(Peak_list peak_list);

extern void unselect_selected_peak_list(Peak_list peak_list);

extern CcpnStatus find_peak_list(Peak_list peak_list, int *first, int *last,
		Block_file block_file, Bool have_low, Bool have_high,
		float low, float high, int *buffer, Bool nonadjacent,
		float drop_factor, float *min_linewidth,
                int ndiagonal_exclusions, Diagonal_exclusion diagonal_exclusions,
                int nexcluded_regions, float ***excluded_regions,
		Bool *dim_checked, int ignore_peak, CcpnString error_msg);

extern CcpnStatus search_region_peak_list(Peak_list peak_list,
		float *first, float *last, Bool *allow_aliasing,
		int *npeaks, int **index_list, CcpnString error_msg);

extern void search_nearest_peak_list(Peak_list peak_list,
		int xdim, int ydim, float xscale, float yscale,
		float *first, float *last,
		Bool *allow_aliasing, int *index_peak, float *d2_min);

/* some helper functions */

extern float determine_peak_scale(Peak peak,
		float intensity_max, float volume_max);

extern void determine_peak_list_max(Peak_list peak_list,
		float *intensity_max, float *volume_max);

#endif /* _incl_peak_list */

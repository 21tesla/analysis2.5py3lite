
/*
======================COPYRIGHT/LICENSE START==========================

peak.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_peak
#define _incl_peak

#include "block_file.h"
#include "drawing_funcs.h"

#define GAUSSIAN_METHOD  0
#define LORENTZIAN_METHOD  1
#define NMETHODS  2

typedef struct Peak
{
    int ndim;
    float *position; /* in points, and counts from 1 */
    int *num_aliasing;
    CcpnString text;
    Bool isSelected;
    float intensity;
    float volume;
    float *text_offset; /* from default position, in points */
    float value_offset; /* from default position, for 1D windows */
/*
    float quality;
*/
    float *line_width; /* in points */
    /*float *box_width;*/ /* in points */
}   *Peak;

extern Peak new_peak(int ndim);

extern void delete_peak(Peak peak);

extern Bool get_is_selected_peak(Peak peak);

extern void set_is_selected_peak(Peak peak, Bool isSelected);

/* peak takes copy of text */
extern CcpnStatus set_text_peak(Peak peak, CcpnString text, CcpnString error_msg);

/* peak takes copy of position */
extern void set_position_peak(Peak peak, float *position);

extern void set_text_offset_peak(Peak peak, int dim, float text_offset);

extern void reset_text_offset_peak(Peak peak, int dim);

/* peak takes copy of num_aliasing */
extern void set_num_aliasing_peak(Peak peak, int *num_aliasing);

extern void draw_peak(Peak peak, int xdim, int ydim, float *first, float *last,
		float xpix, float ypix, float xscale, float yscale,
		int symbol, Bool isSymbolDrawn, Bool isTextDrawn,
		Bool isTextPointerDrawn, Bool hasValueAxis,
		int *tile, Drawing_funcs *drawing_funcs, Generic_ptr data);

extern CcpnStatus get_point_value(float *v, int *point,
			Block_file block_file, CcpnString error_msg);

extern CcpnStatus get_intensity_peak(Peak peak, float *intensity,
			Block_file block_file, CcpnString error_msg);

extern void set_intensity_peak(Peak peak, float intensity);

extern CcpnStatus fit_volume_peak(Peak peak, float *volume, Bool *valid,
		int method, Block_file block_file, Bool *dim_done,
		CcpnString error_msg);

extern void set_volume_peak(Peak peak, float volume);

extern CcpnStatus fit_center_peak(Peak peak, float *center, Bool *valid,
		int method, Block_file block_file, Bool *dim_done,
		CcpnString error_msg);

extern Bool is_in_region_peak(Peak peak, float *first, float *last,
		int *npoints, Bool *allow_aliasing, float *d2_array);

extern void find_scaled_region_peak(Peak peak, int xdim, int ydim,
		float scale, float yscale, float *first, float *last);

/* this fits linewidth in points */
extern CcpnStatus fit_linewidth_peak(Peak peak, float *linewidth,
		Block_file block_file, Bool *dim_done, CcpnString error_msg);

/*
extern void set_quality_peak(Peak peak, float quality);
*/

extern CcpnStatus fit_peaks_in_region(int npeaks, Peak *peaks, int method,
                int *first, int *last, Block_file block_file, Bool *dim_done,
                float *params, float *params_dev, CcpnString error_msg);

extern void set_line_width_peak(Peak peak, int dim, float line_width);

/*
extern void set_box_width_peak(Peak peak, int dim, float box_width);
*/

#endif /* _incl_peak */

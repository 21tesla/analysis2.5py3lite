
/*
======================COPYRIGHT/LICENSE START==========================

peak_cluster.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_peak_cluster
#define _incl_peak_cluster

#include "peak.h"
#include "drawing_funcs.h"

typedef enum { PEAK_CLUSTER_MULTIPLET, PEAK_CLUSTER_SHIFT, PEAK_CLUSTER_STATIC, PEAK_CLUSTER_SYMMETRY } PeakClusterType;

typedef struct Peak_cluster
{
    int ndim;
    PeakClusterType cluster_type;
    CcpnString text;
    CcpnString *dim_text;
    int npeaks;
    int nalloc;
    Peak *peaks;
}   *Peak_cluster;

extern Peak_cluster new_peak_cluster(int ndim, int cluster_type);

extern void delete_peak_cluster(Peak_cluster peak_cluster);

extern CcpnStatus add_peak_peak_cluster(Peak_cluster peak_cluster, Peak peak, CcpnString error_msg);

extern void remove_peak_peak_cluster(Peak_cluster peak_cluster, Peak peak);

/* peak_cluster takes ownership of peaks */
extern void set_peaks_peak_cluster(Peak_cluster peak_cluster, int npeaks, int nalloc, Peak *peaks);
 
extern void clear_peaks_peak_cluster(Peak_cluster peak_cluster);
 
/* peak_cluster takes copy of text */
extern CcpnStatus set_text_peak_cluster(Peak_cluster peak_cluster, CcpnString text, CcpnString error_msg);

/* peak_cluster takes copy of text */
extern CcpnStatus set_dim_text_peak_cluster(Peak_cluster peak_cluster, int dim, CcpnString text, CcpnString error_msg);

extern void draw_peak_cluster(Peak_cluster peak_cluster, int xdim, int ydim,
		Drawing_funcs *drawing_funcs, Generic_ptr data);

#endif /* _incl_peak_cluster */


/*
======================COPYRIGHT/LICENSE START==========================

method.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_method
#define _incl_method

#include "macros.h"
#include "types.h"

#define  FIT_VOLUME_GAUSSIAN3_METHOD	0
#define  NFIT_VOLUME_METHODS		1

#define  FIT_CENTER_PARABOLIC_METHOD	0
#define  FIT_CENTER_GAUSSIAN3_METHOD	1
#define  NFIT_CENTER_METHODS		2

extern float fit_volume_gaussian3_method(int ndim,
			 float y, float *ym, float *yp, Bool *dim_done);

extern void fit_center_parabolic_method(int ndim, float *center,
			 float y, float *ym, float *yp, Bool *dim_done);

extern void fit_center_gaussian3_method(int ndim, float *center,
			 float y, float *ym, float *yp, Bool *dim_done);

#endif /* _incl_method */

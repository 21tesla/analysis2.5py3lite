
/*
======================COPYRIGHT/LICENSE START==========================

geometry.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_geometry
#define _incl_geometry

#include "macros.h"
#include "types.h"

/* return |v| */
extern float vector_length(float *v, int n);

/* return v1 . v2 */
extern float inner_product(float *v1, float *v2, int n);

/* set v1 = v2 ^ v3 */
extern void cross_product(float *v1, float *v2, float *v3);

/* normalise v to have unit length */
extern void normalise_vector(float *v, int n);

/* return angle (radians) between v1 and v2 */
extern float vectors_angle(float *v1, float *v2, int n);

/* set matrix = rotation matrix for angle about axis */
extern void rotation_matrix(float **matrix, float *axis, float angle);

/* set matrix = rotation matrix for rotating v1 to point in v2 direction */
extern void rotation_matrix_vector_to_vector(float **matrix, float *v1, float *v2);

#endif /* _incl_geometry */

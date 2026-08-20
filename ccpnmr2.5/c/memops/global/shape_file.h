
/*
======================COPYRIGHT/LICENSE START==========================

shape_file.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_shape_file
#define _incl_shape_file

#include "macros.h"
#include "types.h"

#include "consts.h"

typedef struct Shape_shape
{
    int size;
    int offset;
    float *values;
}   Shape_shape;

typedef struct Shape_component
{
    float amplitude;
    int nshapes;
    Shape_shape *shapes;
}   Shape_component;

typedef struct Shape_file
{
    int ndim;
    int ncomponents;
    int points[MAX_NDIM];		/* size of data */
    Shape_component *components;
    Integer npoints;			/* total number of points */
    Integer cum_points[MAX_NDIM];       /* cumulative points [size] */
}   *Shape_file;

extern Shape_file new_shape_file(int ndim, int ncomponents, int *points);

extern void delete_shape_file(Shape_file shape_file);

extern CcpnStatus set_amplitude_shape_component(Shape_file shape_file,
	int comp, float amplitude, CcpnString error_msg);

extern CcpnStatus set_data_shape_shape(Shape_file shape_file, int comp, int shape,
	int size, int offset, float *values, CcpnString error_msg);

/* get value at point for specified components */
/* if components is NULL then use all components */
extern float get_point_shape_file(Shape_file shape_file, int *point,
	int ncomponents, int *components);

#endif /* _incl_shape_file */

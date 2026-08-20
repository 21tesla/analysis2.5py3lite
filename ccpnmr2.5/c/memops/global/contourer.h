
/*
======================COPYRIGHT/LICENSE START==========================

contourer.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_contourer
#define _incl_contourer

#include "macros.h"
#include "types.h"

typedef struct Contour_vertex
{
    float x[2];
    struct Contour_vertex *v1; /* previous vertex (NULL if none) */
    struct Contour_vertex *v2; /* next vertex (NULL if none) */
    Bool visited; /* = CCPN_FALSE if not visited, CCPN_TRUE if visited */
}   *Contour_vertex;

typedef struct Contour_vertices
{
    int nvertices; /* number of vertices */
    int nalloc; /* size of each block of vertices allocated */
    Contour_vertex *vertex_store; /* vertex store, ceil(nvertices/nalloc) long */
	/* used to speed up allocation of vertices */
        /* allocate lots in one go and then use one after the other */
        /* so access via vertex i = vertex_store[i/nalloc] + (i%nalloc) */
}   *Contour_vertices;

typedef struct Contours
{
    int n; /* the number of levels */
    Contour_vertices *vertices; /* the contour vertices, one for each level */
}   *Contours;

typedef CcpnStatus (*Contour_get_row_func)(void *user_data, float **data,
							CcpnString error_msg);

/* function keeps hold of pointers, it does not copy them */
/* this makes it easier to change values in between calls to calculate_contours */
extern void *new_contourer_info(void *user_data,
	int nlevels, float *levels, int *npoints,
	float *offset, float *scale,
	Contour_get_row_func contour_get_row_func,
	CcpnString error_msg);

extern void delete_contourer_info(void *contourer_info);

/* information pointed to by contourer_info (e.g. offset and npoints)
   can change in between calls to calculate_contours */
extern Contours calculate_contours(void *contourer_info, CcpnString error_msg);

extern void delete_contours(Contours contours);

typedef CcpnStatus (*Contour_chain_func)(void *user_data, int nvertices,
			Contour_vertex first_vertex, CcpnString error_msg);

extern CcpnStatus process_chains(Contour_vertices contour_vertices,
	void *user_data, Contour_chain_func chain_func, CcpnString error_msg);

#endif /* _incl_contourer */

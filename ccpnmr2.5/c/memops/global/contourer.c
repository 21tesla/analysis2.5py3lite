/*
======================COPYRIGHT/LICENSE START==========================

contourer.c: Part of the CcpNmr Analysis program

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
#include "contourer.h"

/* re-write of contour.c to be thread-safe and more in analysis style */

#define  CONTOUR_NALLOC  50  /* allocate vertices in this size bunch */

typedef struct Contourer_info
{
    void *user_data;
    int nlevels;
    float *levels;
    int *npoints;
    float *offset;
    float *scale;
    Contour_get_row_func get_row_func;
    int nalloc;
    Contour_vertex *v_rows; /* for storing pointers to vertices for a given row */
        /* indexed by level and column, so given index can be NULL or not */
	/* can be used over and over */
}   *Contourer_info;

void *new_contourer_info(void *user_data,
	int nlevels, float *levels, int *npoints,
	float *offset, float *scale,
	Contour_get_row_func get_row_func,
	CcpnString error_msg)
{
    Contourer_info contourer_info;

    sprintf(error_msg, "allocating contourer_info");
    MALLOC_NEW(contourer_info, struct Contourer_info, 1);

    contourer_info->user_data = user_data;
    contourer_info->nlevels = nlevels;
    contourer_info->levels = levels;
    contourer_info->npoints = npoints;
    contourer_info->offset = offset;
    contourer_info->scale = scale;
    contourer_info->get_row_func = get_row_func;

    contourer_info->nalloc = 0;
    contourer_info->v_rows = NULL;

    return contourer_info;
}

void delete_contourer_info(void *contourer_info)
{
    Contourer_info c = (Contourer_info) contourer_info;

    if (c)
	FREE(c->v_rows, Contour_vertex);

    FREE(c, struct Contourer_info);
}

static Contour_vertices new_contour_vertices(void)
{
    Contour_vertices contour_vertices;

    MALLOC_NEW(contour_vertices, struct Contour_vertices, 1);

    contour_vertices->nvertices = 0;
    contour_vertices->nalloc = CONTOUR_NALLOC;
    contour_vertices->vertex_store = NULL;

    return contour_vertices;
}

static Contours new_contours(int n, CcpnString error_msg)
{
    int i;
    Contours contours;

    sprintf(error_msg, "allocating contours");
    MALLOC_NEW(contours, struct Contours, 1);
    MALLOC_NEW(contours->vertices, Contour_vertices, n);

    contours->n = n;
    for (i = 0; i < n; i++)
	contours->vertices[i] = NULL;

    for (i = 0; i < n; i++)
    {
	if (!(contours->vertices[i] = new_contour_vertices()))
	{
	    delete_contours(contours);
	    sprintf(error_msg, "allocating contour_vertices");
	    return NULL;
	}
    }

    return contours;
}

static CcpnStatus init_contours(Contourer_info contourer_info, CcpnString error_msg)
{
    int m = contourer_info->nlevels * (contourer_info->npoints[0] - 1);
    int nalloc = contourer_info->nalloc;

    if (m > nalloc)
    {
	FREE(contourer_info->v_rows, Contour_vertex);
	contourer_info->nalloc = 0;

	MALLOC(contourer_info->v_rows, Contour_vertex, m);

	contourer_info->nalloc = m;
    }

    return CCPN_OK;
}

static Contour_vertex new_vertex(Contour_vertices contour_vertices, CcpnString error_msg)
{
    int nvertices = contour_vertices->nvertices;
    int nalloc = contour_vertices->nalloc;
    int nblocks;
    Contour_vertex v;

    nblocks = nvertices / nalloc;

    if (!(nvertices % nalloc)) /* time for a new block of storage */
    {
	sprintf(error_msg, "allocating contour vertex_store");
	if (nblocks == 0)
	{
	    MALLOC_NEW(contour_vertices->vertex_store, Contour_vertex, 1);
	}
	else
	{
	    REALLOC_NEW(contour_vertices->vertex_store, Contour_vertex, nblocks+1);
	}

	MALLOC_NEW(contour_vertices->vertex_store[nblocks], struct Contour_vertex, nalloc);
    }

    v = contour_vertices->vertex_store[nblocks] + (nvertices % nalloc);
    contour_vertices->nvertices++;

    v->v1 = v->v2 = NULL;
    v->visited = CCPN_FALSE;

    return v;
}

#define  INTERPOLATE(a, b)  ((level - (a)) / ((b) - (a)))

static Contour_vertex new_vertex0(Contour_vertices contour_vertices,
				float *offset, float *scale, float level,
				float d1, float d2, int x, int y, CcpnString error_msg)
{
    Contour_vertex v;

    CHECK_NULL(v = new_vertex(contour_vertices, error_msg));

    v->x[0] = offset[0] + scale[0] * (x + INTERPOLATE(d1, d2));
    v->x[1] = offset[1] + scale[1] * y;

    return v;
}

static Contour_vertex new_vertex1(Contour_vertices contour_vertices,
				float *offset, float *scale, float level,
				float d1, float d2, int x, int y, CcpnString error_msg)
{
    Contour_vertex v;

    CHECK_NULL(v = new_vertex(contour_vertices, error_msg));

    v->x[0] = offset[0] + scale[0] * x;
    v->x[1] = offset[1] + scale[1] * (y + INTERPOLATE(d1, d2));

    return v;
}

#define  NEW_VERTEX0(v, d1, d2, x, y) \
         if (!(v = new_vertex0(contour_vertices, offset, scale, level, d1, d2, x, y, error_msg))) \
	     return CCPN_ERROR;

#define  NEW_VERTEX1(v, d1, d2, x, y) \
         if (!(v = new_vertex1(contour_vertices, offset, scale, level, d1, d2, x, y, error_msg))) \
	     return CCPN_ERROR;

typedef CcpnStatus (*New_edge_func)(Contour_vertices contour_vertices,
		float *offset, float *scale,
		float level, float *data_old, float *data_new,
		Contour_vertex *v_row, Contour_vertex *p_v_col,
		int x, int y, CcpnString error_msg);

/* 0 0
   0 1 */
static CcpnStatus new_edge01(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_old;

    NEW_VERTEX1(v_col, *data_old, *data_new, x+1, y);

    v_old = v_row[x];
    v_old->v1 = v_col;
    v_col->v2 = v_old;
    *p_v_col = v_col;

    return  CCPN_OK;
}

/* 1 1
   1 0 */
static CcpnStatus new_edge32(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_old;

    NEW_VERTEX1(v_col, *data_old, *data_new, x+1, y);

    v_old = v_row[x];
    v_old->v2 = v_col;
    v_col->v1 = v_old;
    *p_v_col = v_col;

    return  CCPN_OK;
}

/* 0 1
   0 0 */
static CcpnStatus new_edge02(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);
    NEW_VERTEX1(v_col, *data_old, *data_new, x+1, y);

    v_row[x] = v_new;
    v_new->v2 = v_col;
    v_col->v1 = v_new;
    *p_v_col = v_col;

    return  CCPN_OK;
}

/* 1 0
   1 1 */
static CcpnStatus new_edge31(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);
    NEW_VERTEX1(v_col, *data_old, *data_new, x+1, y);

    v_row[x] = v_new;
    v_new->v1 = v_col;
    v_col->v2 = v_new;
    *p_v_col = v_col;

    return  CCPN_OK;
}

/* 0 1
   0 1 */
static CcpnStatus new_edge03(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_old, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);

    v_old = v_row[x];
    v_row[x] = v_new;
    v_old->v1 = v_new;
    v_new->v2 = v_old;

    return  CCPN_OK;
}

/* 1 0
   1 0 */
static CcpnStatus new_edge30(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_old, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);

    v_old = v_row[x];
    v_row[x] = v_new;
    v_old->v2 = v_new;
    v_new->v1 = v_old;

    return  CCPN_OK;
}

/* 0 0
   1 0 */
static CcpnStatus new_edge10(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_old;

    v_col = *p_v_col;
    v_old = v_row[x];
    v_old->v2 = v_col;
    v_col->v1 = v_old;

    return  CCPN_OK;
}

/* 1 1
   0 1 */
static CcpnStatus new_edge23(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_old;

    v_col = *p_v_col;
    v_old = v_row[x];
    v_old->v1 = v_col;
    v_col->v2 = v_old;

    return  CCPN_OK;
}

/* 0 0
   1 1 */
static CcpnStatus new_edge11(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX1(v_new, *data_old, *data_new, x+1, y);

    v_col = *p_v_col;
    v_col->v1 = v_new;
    v_new->v2 = v_col;
    *p_v_col = v_new;

    return  CCPN_OK;
}

/* 1 1
   0 0 */
static CcpnStatus new_edge22(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX1(v_new, *data_old, *data_new, x+1, y);

    v_col = *p_v_col;
    v_col->v2 = v_new;
    v_new->v1 = v_col;
    *p_v_col = v_new;

    return  CCPN_OK;
}

/* 0 1
   1 0 */
static CcpnStatus new_edge12(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new, v, v_old;
    float d, d1, d2, d3, d4;

    d1 = *(data_old-1);
    d2 = *data_old;
    d3 = *(data_new-1);
    d4 = *data_new;

    NEW_VERTEX0(v, d3, d4, x, y+1);
    NEW_VERTEX1(v_new, d2, d4, x+1, y);

    d = (d1 + d2 + d3 + d4) / 4;

    v_col = *p_v_col;
    v_old = v_row[x];
    if (d > level)
    {
    	v_col->v1 = v;
    	v->v2 = v_col;
	v_new->v1 = v_old;
	v_old->v2 = v_new;
    }
    else
    {
    	v_col->v1 = v_old;
	v_old->v2 = v_col;
	v_new->v1 = v;
    	v->v2 = v_new;
    }

    v_row[x] = v;
    *p_v_col = v_new;

    return  CCPN_OK;
}

/* 1 0
   0 1 */
static CcpnStatus new_edge21(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new, v, v_old;
    float d, d1, d2, d3, d4;

    d1 = *(data_old-1);
    d2 = *data_old;
    d3 = *(data_new-1);
    d4 = *data_new;

    NEW_VERTEX0(v, d3, d4, x, y+1);
    NEW_VERTEX1(v_new, d2, d4, x+1, y);

    d = (d1 + d2 + d3 + d4) / 4;

    v_col = *p_v_col;
    v_old = v_row[x];
    if (d > level)
    {
    	v_col->v2 = v_old;
	v_old->v1 = v_col;
	v_new->v2 = v;
    	v->v1 = v_new;
    }
    else
    {
    	v_col->v2 = v;
    	v->v1 = v_col;
	v_new->v2 = v_old;
	v_old->v1 = v_new;
    }

    v_row[x] = v;
    *p_v_col = v_new;

    return  CCPN_OK;
}

/* 0 1
   1 1 */
static CcpnStatus new_edge13(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);

    v_col = *p_v_col;
    v_row[x] = v_new;
    v_col->v1 = v_new;
    v_new->v2 = v_col;

    return  CCPN_OK;
}

/* 1 0
   0 0 */
static CcpnStatus new_edge20(Contour_vertices contour_vertices, float *offset, float *scale,
		float level, float *data_old, float *data_new, Contour_vertex *v_row,
		Contour_vertex *p_v_col, int x, int y, CcpnString error_msg)
{
    Contour_vertex v_col, v_new;

    NEW_VERTEX0(v_new, *(data_new-1), *data_new, x, y+1);

    v_col = *p_v_col;
    v_row[x] = v_new;
    v_col->v2 = v_new;
    v_new->v1 = v_col;

    return  CCPN_OK;
}

#define  N  4

#define  GET_DATA_ROW(data) \
	 CHECK_STATUS((get_row_func)(user_data, &data, error_msg));

#define  DATA_ABOVE_LEVEL(d)   ( ((d) > level) ? 1 : 0 )
#define  DATA_ABOVE_LEVEL2(d)  ( ((d) > level) ? 2 : 0 )

static CcpnStatus do_contours(Contourer_info contourer_info,
				Contours contours, CcpnString error_msg)
{
    int *npoints = contourer_info->npoints;
    int nlevels = contourer_info->nlevels;
    float *levels = contourer_info->levels;
    float *offset = contourer_info->offset;
    float *scale = contourer_info->scale;
    Contour_vertex *v_rows;
    void *user_data = contourer_info->user_data;
    Contour_get_row_func get_row_func = contourer_info->get_row_func;
    int l, i0, i1;
    Bool b_old, b_new;
    float level;
    float *data_old, *data_new, *d_old, *d_new, *d_end;
    Contour_vertex *v_row, v_col;
    Contour_vertices contour_vertices;
    static New_edge_func *new_edge;
    static New_edge_func new_edge_func[N][N] =
        { { NULL,       new_edge01, new_edge02, new_edge03 },
          { new_edge10, new_edge11, new_edge12, new_edge13 },
          { new_edge20, new_edge21, new_edge22, new_edge23 },
          { new_edge30, new_edge31, new_edge32, NULL } };

    /* if npoints < 2 then no contours, but get_row_func may need to be called */
    /* so only return early if npoints < 1 in x or y */
    if ((contourer_info->npoints[0] < 1) || (contourer_info->npoints[1] < 1))
	return CCPN_OK;

    CHECK_STATUS(init_contours(contourer_info, error_msg));
    v_rows = contourer_info->v_rows;

    /* first do vertices along bottom row */

    GET_DATA_ROW(data_old);

    i1 = 0;
    for (l = 0; l < nlevels; l++)
    {
	level = levels[l];
	contour_vertices = contours->vertices[l];
	v_row = v_rows + l*(npoints[0]-1);

        d_old = data_old;
        b_old = DATA_ABOVE_LEVEL(*d_old);

        for (i0 = 0; i0 < npoints[0]-1; i0++)
        {
            b_new = DATA_ABOVE_LEVEL(*(++d_old));

            if (b_old ^ b_new) /* i.e. b_old != b_new */
            {
                NEW_VERTEX0(v_row[i0], *(d_old-1), *d_old, i0, i1);
                b_old = b_new;
            }
        }

    }

    for (i1 = 0; i1 < npoints[1]-1; i1++)
    {
        GET_DATA_ROW(data_new);

        for (l = 0; l < nlevels; l++)
        {
            level = levels[l];
	    contour_vertices = contours->vertices[l];
            v_row = v_rows + l*(npoints[0]-1);

            d_old = data_old;
            d_new = data_new;

            d_end = d_old + (npoints[0] - 1);

            b_old = DATA_ABOVE_LEVEL(*d_old)
                        + DATA_ABOVE_LEVEL2(*d_new);

            i0 = 0;
            if ((b_old == 1) || (b_old == 2))
                NEW_VERTEX1(v_col, *d_old, *d_new, i0, i1);

            new_edge = new_edge_func[b_old];
            while (d_old < d_end)
            {
                b_new = DATA_ABOVE_LEVEL(*(++d_old))
                                | DATA_ABOVE_LEVEL2(*(++d_new));

                if (new_edge[b_new])
                {
                    i0 = d_old - data_old - 1;
                    CHECK_STATUS((new_edge[b_new])(contour_vertices, offset, scale, level, d_old, d_new, v_row, &v_col, i0, i1, error_msg));
                    b_old = b_new;
                    new_edge = new_edge_func[b_old];
                }
            }
        }

        SWAP(data_old, data_new, float *);
    }

    return CCPN_OK;
}

Contours calculate_contours(void *contourer_info, CcpnString error_msg)
{
    Contourer_info c = (Contourer_info) contourer_info;
    Contours contours;
    CcpnStatus status;

    CHECK_NULL(contours = new_contours(c->nlevels, error_msg));

    status = do_contours(c, contours, error_msg);

    if (status == CCPN_ERROR)
    {
	delete_contours(contours);
	contours = NULL;
    }

    return contours;
}

static void delete_contour_vertices(Contour_vertices vertices)
{
    int i, nblocks;

    if (!vertices)
	return;

    if (vertices->nvertices > 0)
    {
	nblocks = 1 + (vertices->nvertices-1) / vertices->nalloc;

	for (i = 0; i < nblocks; i++)
	    FREE(vertices->vertex_store[i], struct Contour_vertex);

	FREE(vertices->vertex_store, Contour_vertex);
    }

    FREE(vertices, struct Contour_vertices);
}

void delete_contours(Contours contours)
{
    int i;

    if (!contours)
	return;

    for (i = 0; i < contours->n; i++)
	delete_contour_vertices(contours->vertices[i]);

    FREE(contours->vertices, Contour_vertices);
    FREE(contours, struct Contours);
}

static CcpnStatus process_chain(Contour_vertex v, void *user_data,
			Contour_chain_func chain_func, CcpnString error_msg)
{
    int nvertices;
    Contour_vertex vv;

    nvertices = 1;
    for (vv = v; vv->v1 && (vv->v1 != v); vv = vv->v1)
    {
	nvertices++;
	vv->visited = CCPN_TRUE;
    }

    vv->visited = CCPN_TRUE;

    for (v = v->v2; v && (v != vv); v = v->v2)
    {
	nvertices++;
	v->visited = CCPN_TRUE;
    }

    CHECK_STATUS((chain_func)(user_data, nvertices, vv, error_msg));

    return CCPN_OK;
}

CcpnStatus process_chains(Contour_vertices contour_vertices, void *user_data,
			Contour_chain_func chain_func, CcpnString error_msg)
{
    int i;
    int nvertices = contour_vertices->nvertices;
    int nalloc = contour_vertices->nalloc;
    Contour_vertex v;
    Contour_vertex *vertex_store = contour_vertices->vertex_store;

    for (i = 0; i < nvertices; i++)
    {
	v = vertex_store[i/nalloc] + i % nalloc;
	v->visited = CCPN_FALSE;
    }

    for (i = 0; i < nvertices; i++)
    {
	v = vertex_store[i/nalloc] + i % nalloc;

	if (v->visited)
	    continue;

	CHECK_STATUS(process_chain(v, user_data, chain_func, error_msg));
    }

    return CCPN_OK;
}

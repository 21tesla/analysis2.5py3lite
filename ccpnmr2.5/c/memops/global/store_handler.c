
/*
======================COPYRIGHT/LICENSE START==========================

store_handler.c: Part of the CcpNmr Analysis program

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
#include "store_handler.h"

#include "store_consts.h"
#include "utility.h"

#define  CHECK_FOR_SWAP(x, n) \
	 if (store_handler->swap) \
	     swap_bytes((char *) (x), BYTES_PER_WORD*n);

#define  STORE_FWRITE(x, n) \
	 {   CHECK_FOR_SWAP((x), n); \
	     fwrite((void *) (x), BYTES_PER_WORD, n, store_handler->fp); \
	     CHECK_FOR_SWAP((x), n);   }

static void store_start_draw(Generic_ptr data)
{
}

static void store_end_draw(Generic_ptr data)
{
}

static void store_new_draw_range(Generic_ptr data, float x0, float y0,
				float x1, float y1, Bool clip)
{
}

static void store_draw_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
}

static void store_draw_clipped_line(Generic_ptr data, float x0, float y0, float x1, float y1)
{
}

static void store_draw_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_polyline_store_handler((Store_handler) data, polyline);
}

static void store_draw_clipped_polyline(Generic_ptr data, Poly_line polyline)
{
    draw_clipped_polyline_store_handler((Store_handler) data, polyline);
}

static void store_draw_text(Generic_ptr data, CcpnString text, float x, float y,
                                                        float a, float b)
{
}
 
static void store_set_draw_color(Generic_ptr data, float *color)
{
}
 
static void store_set_draw_font(Generic_ptr data, CcpnString name, int size)
{
}
 
static void store_set_line_style(Generic_ptr data, int line_style)
{
}
 
static void store_set_line_width(Generic_ptr data, float line_width)
{
}

static void store_fill_circle(Generic_ptr data, float x, float y, float r)
{
}

static void store_draw_circle(Generic_ptr data, float x, float y, float r)
{
}

static void store_get_background(Generic_ptr data, float *color)
{
}

static void store_get_region(Generic_ptr data, float *x0, float *y0,
                                                    float *x1, float *y1)
{
}

static Drawing_funcs drawing_funcs = { STORE_DRAWING,
	store_start_draw, store_end_draw, store_new_draw_range,
	store_draw_line, store_draw_clipped_line,
	store_draw_polyline, store_draw_clipped_polyline, store_draw_text,
	store_set_draw_color, store_set_draw_font, store_set_line_style,
        store_set_line_width, store_fill_circle, store_draw_circle,
        store_get_background, store_get_region };
 
Drawing_funcs *store_drawing_funcs()
{
    return &drawing_funcs;
}

static int get_dir_index(Store_handler store_handler)
{
    int array[MAX_NDIM];
    int a, i, dir_index, ndim = store_handler->ndim;

    for (i = 0; i < ndim; i++)
    {
	if ((i == store_handler->xdim) || (i == store_handler->ydim))
	    a = store_handler->block[i];
	else if (store_handler->block[i] == 0)
	    a = store_handler->plane[i];
	else
	    a = store_handler->block_size[i] * store_handler->block[i] +
			store_handler->plane[i] - store_handler->first_plane[i];

	array[i] = a;
    }

    INDEX_OF_ARRAY(dir_index, array, store_handler->cumul_dir, ndim);

    if (store_handler->have_pos && store_handler->have_neg)
    {
	dir_index *= 2;
/* relies on negative contours being before positive ones */
	if (store_handler->is_pos)
	    dir_index++;
    }

    return dir_index;
}

Store_handler new_store_handler(CcpnString file_name,
					Bool swap, CcpnString error_msg)
{
    FILE *fp;
    Store_handler store_handler;

    sprintf(error_msg, "allocating store_handler");
    MALLOC_NEW(store_handler, struct Store_handler, 1);

    if (OPEN_FOR_BINARY_WRITING(fp, file_name))
    {
	sprintf(error_msg, "could not open '%s' for writing", file_name);
	FREE(store_handler, struct Store_handler);
	return NULL;
    }

    store_handler->fp = fp;
    store_handler->swap = swap;
    store_handler->directory = NULL;

    return store_handler;
}

void delete_store_handler(Store_handler store_handler)
{
    FCLOSE(store_handler->fp);

    FREE(store_handler, struct Store_handler);
}

void draw_polyline_store_handler(Store_handler store_handler, Poly_line polyline)
{
    int i, ind, m, n = polyline->nvertices, npolys = 0;
    Point2f *v = polyline->vertices;
    int dir_index = get_dir_index(store_handler);

    if (store_handler->directory[dir_index] == -1)
    {
	store_handler->directory[dir_index] = store_handler->offset;
	STORE_FWRITE(&npolys, 1);  /* place holder */
        store_handler->offset++;
    }

    if (polyline->closed)
	m = -n;
    else
	m = n;

    STORE_FWRITE(&m, 1);

    for (i = 0; i < n; i++)
    {
	STORE_FWRITE(&(v[i].x), 1);
	STORE_FWRITE(&(v[i].y), 1);
    }

    store_handler->offset += 1 + 2*n;
    if (store_handler->is_pos)
	store_handler->npos_polys++;
    else
	store_handler->nneg_polys++;
}

void draw_clipped_polyline_store_handler(Store_handler store_handler, Poly_line polyline)
{
    draw_polyline_store_handler(store_handler, polyline);
}

CcpnStatus init_store_save(Store_handler store_handler,
			int ndim, int xdim, int ydim,
			int *npoints, int *first, int *last,
			int *block_size, int *nblocks,
			int nlevels, float *levels,
			CcpnString error_msg)
{
    int i, n, p, dir_size, header_size;
    int *directory, *header;
    float *hdr;
    Bool have_pos, have_neg;

/* must be kept in sync with store_file.c */

/* header size not including directories */
    header_size = HEADER0 + ndim*HEADER1 + HEADER2*nlevels;
    sprintf(error_msg, "allocating header");
    MALLOC(header, int, header_size);

    n = 0;
    header[n++] = STORE_MAGIC;  /* magic constant and also determines byte ordering */
    header[n++] = STORE_VERSION;  /* version number */
    header[n++] = header_size;
    header[n++] = ndim;
    header[n++] = xdim;
    header[n++] = ydim;

    COPY_VECTOR(header+n, npoints, ndim);
    n += ndim;

    COPY_VECTOR(header+n, first, ndim);
    n += ndim;

    COPY_VECTOR(header+n, last, ndim);
    n += ndim;

    COPY_VECTOR(header+n, block_size, ndim);
    n += ndim;

    header[n++] = nlevels;
    hdr = (float *) (header + n);
    COPY_VECTOR(hdr, levels, nlevels);
    n += nlevels;

    if (n != header_size)
    {
	sprintf(error_msg, "inconsistent header size: %d vs %d", n, header_size);
	return CCPN_ERROR;
    }

    have_pos = CCPN_FALSE;
    have_neg = CCPN_FALSE;
    for (i = 0; i < nlevels; i++)
    {
	if (levels[i] >= 0)
	    have_pos = CCPN_TRUE;
	else
	    have_neg = CCPN_TRUE;
    }

    dir_size = 1;
    for (i = 0; i < ndim; i++)
    {
	if ((i == xdim) || (i == ydim))
	    p = nblocks[i];
	else
	    p = last[i] - first[i];

	store_handler->cumul_dir[i] = dir_size;
	dir_size *= p;
    }

    if (have_pos && have_neg)
    	dir_size *= 2;

/*
    printf("dir_size = %d\n", dir_size);
*/

    sprintf(error_msg, "allocating contour store directory");
    MALLOC(directory, int, dir_size);

    for (i = 0; i < dir_size; i++)
	directory[i] = -1;

    CHECK_FOR_SWAP(header, n);
    if (FWRITE(header, BYTES_PER_WORD, n, store_handler->fp))
	RETURN_ERROR_MSG("writing contour store header");
    CHECK_FOR_SWAP(header, n);

    CHECK_FOR_SWAP(directory, dir_size);
    if (FWRITE(directory, BYTES_PER_WORD, dir_size, store_handler->fp))
	RETURN_ERROR_MSG("writing contour store directory");
    CHECK_FOR_SWAP(directory, dir_size);

    store_handler->ndim = ndim;
    store_handler->xdim = xdim;
    store_handler->ydim = ydim;
    store_handler->header_size = header_size;
    store_handler->offset = 0;
    store_handler->have_pos = have_pos;
    store_handler->have_neg = have_neg;
    store_handler->dir_size = dir_size;
    store_handler->directory = directory;
    COPY_VECTOR(store_handler->block_size, block_size, ndim);

    FREE(header, int);

    return CCPN_OK;
}

/* note that this block is relative to block_min */
void init_store_block(Store_handler store_handler, int *block)
{
    COPY_VECTOR(store_handler->block, block, store_handler->ndim);
}

void init_store_plane(Store_handler store_handler, int *plane)
{
    COPY_VECTOR(store_handler->plane, plane, store_handler->ndim);
    store_handler->npos_polys = 0;
    store_handler->nneg_polys = 0;
}

void init_store_level(Store_handler store_handler, float level)
{
    store_handler->is_pos = (level >= 0 ? CCPN_TRUE : CCPN_FALSE);
}

void end_store_plane(Store_handler store_handler)
{
    int dir_index, offset;

    if (store_handler->nneg_polys > 0)
    {
	store_handler->is_pos = CCPN_FALSE;
	dir_index = get_dir_index(store_handler);
	offset = store_handler->directory[dir_index];
	offset += store_handler->header_size + store_handler->dir_size;
	FSEEK_ABSOLUTE(store_handler->fp, offset);
	STORE_FWRITE(&store_handler->nneg_polys, 1);
    }

    if (store_handler->npos_polys > 0)
    {
	store_handler->is_pos = CCPN_TRUE;
	dir_index = get_dir_index(store_handler);
	offset = store_handler->directory[dir_index];
	offset += store_handler->header_size + store_handler->dir_size;
	FSEEK_ABSOLUTE(store_handler->fp, offset);
	STORE_FWRITE(&store_handler->npos_polys, 1);
    }

    if (store_handler->npos_polys > 0 || store_handler->nneg_polys > 0)
    {
	offset = store_handler->offset;
	offset += store_handler->header_size + store_handler->dir_size;
	FSEEK_ABSOLUTE(store_handler->fp, offset);
    }
}

void end_store_save(Store_handler store_handler)
{
    FSEEK_ABSOLUTE(store_handler->fp, store_handler->header_size);

    STORE_FWRITE(store_handler->directory, store_handler->dir_size);

    FCLOSE(store_handler->fp);

/*
    { int i, ind;
    for (i = 0; i < store_handler->dir_size; i++)
    {
        ind = store_handler->directory[i];
	if (ind != -1)
	    printf("%d: %d\n", i, ind);
    } }
*/
}


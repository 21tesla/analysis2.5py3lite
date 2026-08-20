
/*
======================COPYRIGHT/LICENSE START==========================

store_file.c: Part of the CcpNmr Analysis program

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
#include "store_file.h"

#include "store_consts.h"
#include "utility.h"

#define  CHECK_FOR_SWAP(x, n) \
         if (store_file->swap) \
             swap_bytes((char *) (x), BYTES_PER_WORD*n);

#define  STORE_FREAD(x, n, err_msg) \
	 {   if (FREAD(x, BYTES_PER_WORD, n, store_file->fp)) \
	     {   strcpy(error_msg, err_msg); \
	         return CCPN_ERROR;   } \
	     CHECK_FOR_SWAP(x, n);   }

#define  CHECK_EXPECTED_VALUE(x, y, msg) \
	 {   if ((x) != (y)) \
	     {   sprintf(error_msg, "%s: got %d, expected %d", msg, x, y); \
		 return CCPN_ERROR;   }   }

#define  CHECK_EXPECTED_VALUES(x, y, msg, n) \
	 {   int I; \
	     for (I = 0; I < (n); I++) \
	     {   Line Msg; \
		 sprintf(Msg, "%s (element %d)", msg, I+1); \
		 CHECK_EXPECTED_VALUE((x)[I], (y)[I], Msg); \
	     }   }

static int get_dir_index(Store_file store_file, int *block, int *plane)
{
    int array[MAX_NDIM];
    int a, i, dir_index, ndim = store_file->ndim;

    for (i = 0; i < ndim; i++)
    {
	if ((i == store_file->xdim) || (i == store_file->ydim))
	    a = block[i];
	else
	    a = store_file->block_size[i] * block[i] + plane[i];

	array[i] = a;
    }

    INDEX_OF_ARRAY(dir_index, array, store_file->cumul_dir, ndim);

    if (store_file->have_pos && store_file->have_neg)
	dir_index *= 2;

    return dir_index;
}

static CcpnStatus init_store_file(Store_file store_file, CcpnString file_name,
							CcpnString error_msg)
{
    int t, ndim, nlevels, i, p, bmin, bmax, dir_size;
    int block_size[MAX_NDIM];
    float level;
    FILE *fp;

/* must be kept in sync with store_handler.c */

    CHECK_OPEN_FOR_BINARY_READING(fp, file_name);

    store_file->fp = fp;

    if (FREAD(&t, BYTES_PER_WORD, 1, fp))
	RETURN_ERROR_MSG("reading magic word (file empty?)");

    store_file->swap = CCPN_FALSE;
    if (t != STORE_MAGIC)
    {
	swap_bytes((char *) &t, BYTES_PER_WORD);
	if (t != STORE_MAGIC)
	    RETURN_ERROR_MSG("first word of file is not magic word (normal order or byte swapped)");

	store_file->swap = CCPN_TRUE;
    }

    STORE_FREAD(&t, 1, "reading version");
    CHECK_EXPECTED_VALUE(t, STORE_VERSION, "version number");
/*
    printf("version = %d\n", t);
*/

    STORE_FREAD(&store_file->header_size, 1, "reading header_size");
/*
    printf("header_size = %d\n", store_file->header_size);
*/

    STORE_FREAD(&ndim, 1, "reading ndim");
    CHECK_EXPECTED_VALUE(ndim, store_file->ndim, "ndim");
/*
    printf("ndim = %d\n", ndim);
*/

    STORE_FREAD(&t, 1, "reading xdim");
    CHECK_EXPECTED_VALUE(t, store_file->xdim, "xdim");
/*
    printf("xdim = %d\n", t);
*/

    STORE_FREAD(&t, 1, "reading ydim");
    CHECK_EXPECTED_VALUE(t, store_file->ydim, "ydim");
/*
    printf("ydim = %d\n", t);
*/

    STORE_FREAD(store_file->npoints, ndim, "reading npoints");
    STORE_FREAD(store_file->first, ndim, "reading first");
    STORE_FREAD(store_file->last, ndim, "reading last");

    STORE_FREAD(block_size, ndim, "reading block_size");
    CHECK_EXPECTED_VALUES(block_size, store_file->block_size, "block_size", ndim);

    STORE_FREAD(&nlevels, 1, "reading nlevels");

    for (i = 0; i < nlevels; i++)
    {
        STORE_FREAD(&level, 1, "reading level");
	if (level >= 0)
	    store_file->have_pos = CCPN_TRUE;
	else
	    store_file->have_neg = CCPN_TRUE;
    }

    dir_size = 1;
    for (i = 0; i < ndim; i++)
    {
        bmin = store_file->first[i] / block_size[i];
        if ((i == store_file->xdim) || (i == store_file->ydim))
	{
            bmax = (store_file->last[i]-1) / block_size[i];
            p = bmax - bmin + 1;
	}
        else
	{
            p = store_file->last[i] - store_file->first[i];
	}

        store_file->cumul_dir[i] = dir_size;
        dir_size *= p;

	store_file->block_min[i] = bmin;
	store_file->nblocks[i] = 1 + (store_file->npoints[i] - 1) / store_file->block_size[i];
    }

    if (store_file->have_pos && store_file->have_neg)
        dir_size *= 2;

/*
    printf("dir_size = %d\n", dir_size);
*/

    store_file->dir_size = dir_size;

    sprintf(error_msg, "allocating contour store directory");
    MALLOC(store_file->directory, int, dir_size);
    STORE_FREAD(store_file->directory, dir_size, "reading directory");

/*
    for (i = 0; i < dir_size; i++)
	if (store_file->directory[i] != -1)
	    printf("directory[%d] = %d\n", i, store_file->directory[i]);
*/

    return CCPN_OK;
}

Store_file new_store_file(CcpnString file_name, int ndim,
	int xdim, int ydim, int *block_size, CcpnString error_msg)
{
    Store_file store_file;

    sprintf(error_msg, "allocating store_file");
    MALLOC_NEW(store_file, struct Store_file, 1);

    store_file->fp = NULL;
    store_file->ndim = ndim;
    store_file->xdim = xdim;
    store_file->ydim = ydim;
    COPY_VECTOR(store_file->block_size, block_size, ndim);

    store_file->have_pos = CCPN_FALSE;
    store_file->have_neg = CCPN_FALSE;
    store_file->directory = NULL;

    if (init_store_file(store_file, file_name, error_msg) == CCPN_ERROR)
    {
	delete_store_file(store_file);
	return NULL;
    }

    return store_file;
}

void delete_store_file(Store_file store_file)
{
    FCLOSE(store_file->fp);

    FREE(store_file, struct Store_file);
}

static CcpnStatus read_polylines(Store_file store_file,
			int npolylines, Poly_line *polylines,
			Bool transposed, CcpnString error_msg)
{
    int i, j, m, n;
    Poly_line polyline;
    Point2f *v;

    for (i = 0; i < npolylines; i++)
    {
	polyline = polylines[i];

/* TBD: should put in block, plane, level?? */
	sprintf(error_msg, "reading number of vertices of polyline %d", i+1);
	STORE_FREAD(&m, 1, error_msg);
	n = ABS(m);

	sprintf(error_msg, "allocating vertices memory");
	MALLOC(polyline->vertices, Point2f, n);
	polyline->closed = (m > 0) ? CCPN_FALSE : CCPN_TRUE;
	polyline->nvertices = n;

	for (j = 0; j < n; j++)
	{
	    v = polyline->vertices + j;
	    sprintf(error_msg, "reading polyline %d vertex %d", i+1, j+1);
	    STORE_FREAD(&(v->x), 1, error_msg);
	    STORE_FREAD(&(v->y), 1, error_msg);
	    if (transposed)
		SWAP(v->x, v->y, float);
	}
    }

    return CCPN_OK;
}

static CcpnStatus process_contours(Store_file store_file,
		int *block, int *plane, int dir_index,
		int level, Store_poly_func poly_func,
		void *user_data, Bool transposed, CcpnString error_msg)
{
    int i, n, npolylines;
    int offset = store_file->directory[dir_index];
    FILE *fp = store_file->fp;
    Poly_line *polylines;

    if (offset >= 0)
    {
	offset += store_file->header_size + store_file->dir_size;

	FSEEK_ABSOLUTE(fp, offset);

	STORE_FREAD(&npolylines, 1, "reading number of polylines");

	sprintf(error_msg, "allocating polyline memory");
	MALLOC(polylines, Poly_line, npolylines);
	for (i = 0; i < npolylines; i++)
	    polylines[i] = NULL;

	for (i = 0; i < npolylines; i++)
	{
	    MALLOC(polylines[i], struct Poly_line, 1);
	    polylines[i]->vertices = (Point2f *) NULL;
	}

	if (read_polylines(store_file, npolylines, polylines, transposed, error_msg) == CCPN_ERROR)
	{
	    for (i = 0; i < npolylines; i++)
	    {
	    FREE(polylines[i]->vertices, Point2f);
	    FREE(polylines[i], struct Poly_line);
	    }

	    FREE(polylines, Poly_line);

	    return CCPN_ERROR;
	}
    }
    else
    {
	npolylines = 0;
	polylines = NULL;
    }

    return (poly_func)(user_data, block, plane, level, npolylines, polylines, error_msg);
}

/* note that this block is relative to 0, not block_min */
CcpnStatus process_contours_store_file(Store_file store_file,
		int *block, int *plane, Store_poly_func poly_func,
		void *user_data, Bool transposed, CcpnString error_msg)
{
    int blk[MAX_NDIM], pln[MAX_NDIM];
    int i, p, dir_index, offset, ndim = store_file->ndim;

    for (i = 0; i < ndim; i++)
    {
	blk[i] = block[i] - store_file->block_min[i];

        if ((i != store_file->xdim) && (i != store_file->ydim))
        {
            p = block[i] * store_file->block_size[i];
	    p = store_file->first[i] - p;
	    p = MAX(0, p);
            pln[i] = plane[i] - p;
        }
        else
        {
            pln[i] = 0;
        }
    }

    /* below is the index for negative if have negative */
    /* but is the index for positive if do not have negative */
    dir_index = get_dir_index(store_file, blk, pln);
    CHECK_STATUS(process_contours(store_file, block, plane, dir_index, 0, poly_func, user_data, transposed, error_msg));

    if (store_file->have_pos && store_file->have_neg)
    {
        /* below is the index for positive if have negative */
	dir_index++;
	CHECK_STATUS(process_contours(store_file, block, plane, dir_index, 1, poly_func, user_data, transposed, error_msg));
    }

    return CCPN_OK;
}
 

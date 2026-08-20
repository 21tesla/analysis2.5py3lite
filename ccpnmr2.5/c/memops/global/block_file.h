
/*
======================COPYRIGHT/LICENSE START==========================

block_file.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_block_file
#define _incl_block_file

/* designed for random access read of blocked data files */

#include "macros.h"
#include "types.h"

#include "consts.h"
#include "hash_table.h"
#include "int_array.h"
#include "mem_cache.h"
#include "shape_file.h"

/* all counts are in terms of real numbers even for complex data */

typedef struct Block_data
{
    Int_array block;			/* the block (key so needs to be Int_array rather than Integer*) */
    int ndim;				/* number of dimensions */
    Integer size;			/* size in points (disk and memory) */
    int block_size[MAX_NDIM];		/* size of block by dimension */
    Integer cum_block_size[MAX_NDIM];	/* cumulative size of block by dimension */
    Float *data;
}   *Block_data;

typedef enum {BLOCK_FILE, SHAPE_FILE} Block_file_kind;

#define  BLOCK_FILE_ATTRS \
    Block_file_kind block_file_kind; \
    char *file; \
    FILE *fp; \
    int ndim; \
    int points[MAX_NDIM];		/* size of data */ \
    int block_size[MAX_NDIM];		/* size of block (except possibly for end blocks) */ \
    Bool dim_wrapped[MAX_NDIM];		/* whether the data should be wrapped */ \
    Bool padded;			/* whether end blocks are padded */ \
    Bool writeable;			/* whether file is writeable or only readable */ \
    Integer npoints;			/* total number of points */ \
    Integer nblocks;			/* total number of blocks */ \
    int blocks[MAX_NDIM];		/* number of blocks =~ size / block_size */ \
    Integer cum_points[MAX_NDIM];	/* cumulative points [size] */ \
    Integer cum_blocks[MAX_NDIM];	/* cumulative blocks [blocks] */ \
    Integer cum_block_points[MAX_NDIM];	/* cumulative block points [padded size, block size] */

typedef struct Block_file
{
    BLOCK_FILE_ATTRS
    int bytes_per_point;		/* number of bytes per data point */
    Bool big_endian;			/* whether data is big or little endian */
    int header;				/* length of file header in bytes */
    int block_header;			/* length of each block header in bytes (for Varian) */
    Bool integer;			/* whether data is integer or real */
    Mem_cache mem_cache;
    Hash_table block_table;		/* blocks in cache */
    Hash_table dirty_blocks;		/* blocks that have been modified */
}   *Block_file;

typedef struct Shape_block_file
{
    BLOCK_FILE_ATTRS
    Shape_file shape_file;
}   *Shape_block_file;

/* should not be used for block data associated with a block file */
extern void delete_block_data(Block_data block_data);

/* does not open file */
extern Block_file new_block_file(char *file, int ndim,
	int *points, int *block_size, Bool *dim_wrapped,
	Mem_cache mem_cache, int bytes_per_point, Bool big_endian,
	Bool padded, int header, Bool integer, Bool writeable, int block_header);

/* treat Shape_block_file as a Block_file (sort of a subclass but not quite) */
extern Block_file new_shape_block_file(char *file, int ndim, int *points,
	int *block_size, Bool *dim_wrapped, Shape_file shape_file);

/* closes file if open */
extern void delete_block_file(Block_file block_file);

extern CcpnStatus open_block_file(Block_file block_file);

extern void close_block_file(Block_file block_file);

extern CcpnStatus have_peak_block_file(Block_file block_file,
	int *first, int *last, Bool have_low, Bool have_high,
	float low, float high, Bool *have_peak, CcpnString error_msg);

extern CcpnStatus max_value_block_file(Block_file block_file,
	int *first, int *last, float *max_value, CcpnString error_msg);

extern CcpnStatus min_value_block_file(Block_file block_file,
	int *first, int *last, float *min_value, CcpnString error_msg);

/* this locks returned block_data */
extern Block_data get_block_data_number(Block_file block_file, Integer block);

/* this locks returned block_data */
extern Block_data get_block_data_array(Block_file block_file, int *block);

/* this locks returned block_data if a BLOCK_FILE */
extern Block_data get_components_block_data(Block_file block_file,
	int *block_array, int ncomponents, int *components);

extern CcpnStatus lock_block_data(Block_file block_file, Block_data block_data);

extern CcpnStatus unlock_block_data(Block_file block_file, Block_data block_data);

/* get the value at a (grid) point */
extern CcpnStatus get_point_block_file(Block_file block_file, float *v,
					int *point, CcpnString error_msg);

/* get the value at a (grid) point */
CcpnStatus get_components_point_block_file(Block_file block_file, float *v, int *point,
                        int ncomponents, int *components, CcpnString error_msg);

/* get the value at a (not necessarily grid) point */
extern CcpnStatus get_position_block_file(Block_file block_file, float *v,
				float *position, CcpnString error_msg);

/* Find the linewidth (in points) for a given point in a given dim */
/* find_maximum is CCPN_TRUE if point is a maximum, otherwise it's CCPN_FALSE */
/* Note, v is normally the value of a peak at that point */
/* And normally have_maximum = CCPN_TRUE if v > 0, otherwise it's CCPN_FALSE */
extern CcpnStatus linewidth_block_file(Block_file block_file,
			Bool have_maximum, float v, int *point, int dim,
			float *linewidth, CcpnString error_msg);

/* get the values inside a box, with the values returned as a 1D array */
extern CcpnStatus get_box_block_file(Block_file block_file, int *p_npoints, float **p_v,
                                int *box_min, int *box_max, CcpnString error_msg);

/* set the value at a (grid) point */
extern CcpnStatus set_point_block_file(Block_file block_file, float v,
						int *point, CcpnString error_msg);

/* set the values inside a box, with the values input as a 1D array */
/* does not save the data to disk */
extern CcpnStatus set_box_block_file(Block_file block_file, int npoints, float *v,
                                int *box_min, int *box_max, CcpnString error_msg);

/* saves blocks to disk that have been touched since the last save */
extern CcpnStatus save_block_file(Block_file block_file, CcpnString error_msg);

extern CcpnStatus check_block_file(Block_file block_file, CcpnString error_msg);

/* calculates the normalised inner product between two slices */
/* dim1, dim2 count from 0 but slice1 and slice2 are in points which count from 1 */
/* slice1 and slice2 are of length ndim-1, not ndim (so skipping dim1 and dim2 index) */
/* point 1 in slice1 corresponds with point 1+offset in slice2 */
/* point 2 in slice1 corresponds with point 1+offset+scale in slice2 */
/* etc. */
extern CcpnStatus compare_slices_block_files(Block_file block_file1, int dim1, int *slice1,
		Block_file block_file2, int dim2, int *slice2, float *v,
		float scale, float offset, CcpnString error_msg);

/* xdim = the x dimension for this block_file (counting from 0) */
/* ydim = the y dimension for this block_file (counting from 0) */
/* xmin = lowest grid point in xdim (counting from 0) */
/* xmax = highest grid point in xdim (counting from 0) */
/* ymin = lowest grid point in ydim (counting from 0) */
/* ymax = highest grid point in ydim (counting from 0) */
/* plane = fixed point in orthogonal dimensions */
typedef struct
{
    Block_file block_file;
    int xdim;
    int ydim;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float plane[MAX_NDIM];
}   Block_plane;

extern Block_plane *alloc_block_planes(int nplanes);

extern void free_block_planes(Block_plane *block_planes);

/* calculate the product of nplanes block_files (which are ND, where N can vary) */
/* and put the result in output block_file (which must be 2D) */
extern CcpnStatus plane_product_block_files(Block_file block_file,
		int nplanes, Block_plane *block_planes,
		float scale, CcpnString error_msg);

#endif /* _incl_block_file */


/*
======================COPYRIGHT/LICENSE START==========================

block_file.c: Part of the CcpNmr Analysis program

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
#include "block_file.h"

#include "utility.h"

static Bool equal_block(Hash_key key1, Hash_key key2)
{
    Int_array block1 = (Int_array) key1;
    Int_array block2 = (Int_array) key2;

    return equal_int_array(block1, block2);
}

static Hash_value hash_block(Hash_key key)
{
    Int_array block = (Int_array) key;

    return hash_int_array(block);
}

static void size_block_data(Block_file block_file, Block_data block_data)
{
    int i, s;
    Integer size = 1;

    for (i = 0; i < block_file->ndim; i++)
    {
	s = block_file->points[i] % block_file->block_size[i];
	if (!s || block_file->padded ||
		(block_data->block->values[i] != (block_file->blocks[i]-1)))
	    s = block_file->block_size[i];

	block_data->block_size[i] = s;
	block_data->cum_block_size[i] = size;
	size *= s;
    }

    block_data->size = MIN(size, block_file->npoints);
}

static Block_data new_block_data(Block_file block_file, Int_array block)
{
    Block_data block_data;

    MALLOC_NEW(block_data, struct Block_data, 1);
    block_data->ndim = block_file->ndim;
    block_data->block = block;

    size_block_data(block_file, block_data);
    MALLOC_NEW(block_data->data, Float, block_data->size);

    return block_data;
}

void delete_block_data(Block_data block_data)
{
    delete_int_array(block_data->block);
    FREE(block_data->data, Float);
    FREE(block_data, struct Block_data);
}

static void delete_data_block(Generic_ptr object, Generic_ptr data)
{
    Block_data block_data = (Block_data) object;
    Block_file block_file = (Block_file) data;

    (void) remove_hash_table(block_file->block_table,
					(Hash_key) block_data->block);
    delete_block_data(block_data);
}

static CcpnStatus read_block_data(Block_file block_file, Block_data block_data)
{
    int i, j, k1, k2, m1, m2;
    Integer offset, block;

#ifndef WIN64
    if (!block_file->fp)
	return CCPN_ERROR;
#else
    if (!block_file->fp)
       CHECK_STATUS(open_block_file(block_file));
#endif

    INDEX_OF_ARRAY(offset, block_data->block->values, block_file->cum_block_points, block_file->ndim);
    offset = block_file->bytes_per_point*offset + block_file->header;

    if (block_file->block_header != 0)
    {
	INDEX_OF_ARRAY(block, block_data->block->values, block_file->cum_blocks, block_file->ndim);
	/* the assumption here is that block 0 also needs extra offset
	   so block_file->header does not take that into account, hence +1 */
	offset += (block+1) * block_file->block_header;
    }

/*
    printf("offset = %ld\n", offset);
    for (i = 0; i < block_file->ndim; i++)
	printf(" %d: bs=%d, ps=%d, blk=%d\n", i, block_file->block_size[i], block_file->points[i], block_data->block->values[i]);
*/

/*
    if (FSEEK_ABSOLUTE_BYTES(block_file->fp, offset))
*/
    if (fseek_absolute_bytes(block_file->fp, offset))
	return CCPN_ERROR;

    if (FREAD(block_data->data, block_file->bytes_per_point, block_data->size, block_file->fp))
	return CCPN_ERROR;

    if (block_file->bytes_per_point < BYTES_PER_WORD)
    {
	for (i = block_data->size-1; i >= 0; i--)
	{
	    if (block_file->big_endian)
	    {
		m1 = BYTES_PER_WORD - block_file->bytes_per_point;
		m2 = - block_file->bytes_per_point;
	    }
	    else
	    {
		m1 = m2 = 0;
	    }

	    for (j = 0; j < block_file->bytes_per_point; j++)
	    {
		k1 = i*BYTES_PER_WORD + m1 + j;
		k2 = i*block_file->bytes_per_point + j;
		block_data->data[k1] = block_data->data[k2];
	    }

	    for (j = block_file->bytes_per_point; j < BYTES_PER_WORD; j++)
	    {
		k1 = i*BYTES_PER_WORD + m2 + j;
		block_data->data[k1] = 0;
	    }
	}
    }

    if (block_file->big_endian != is_big_endian())
	swap_bytes((char *) block_data->data, BYTES_PER_WORD*block_data->size);

    if (block_file->integer)
	float_words(block_data->data, block_data->size);

    return CCPN_OK;
}

static CcpnStatus write_block_data(Block_file block_file, Block_data block_data)
{
    Integer offset;
    Bool not_ok;

    if (!block_file->fp)
	return CCPN_ERROR;

    INDEX_OF_ARRAY(offset, block_data->block->values, block_file->cum_block_points, block_file->ndim);
    offset = block_file->bytes_per_point*offset + block_file->header;

/*
    if (FSEEK_ABSOLUTE_BYTES(block_file->fp, offset))
*/
    if (fseek_absolute_bytes(block_file->fp, offset))
	return CCPN_ERROR;

    if (block_file->integer)
	int_words((int *) block_data->data, block_data->size);

    if (block_file->big_endian != is_big_endian())
	swap_bytes((char *) block_data->data, BYTES_PER_WORD*block_data->size);

    not_ok = FWRITE(block_data->data, block_file->bytes_per_point, block_data->size, block_file->fp);

    if (block_file->big_endian != is_big_endian())
	swap_bytes((char *) block_data->data, BYTES_PER_WORD*block_data->size);

    if (block_file->integer)
	float_words(block_data->data, block_data->size);

    if (not_ok)
	return CCPN_ERROR;

    return CCPN_OK;
}

static CcpnStatus init_block_file(Block_file block_file, char *file,
	int ndim, int *points, int *block_size, Bool padded)
{
    int i, s;
    Integer t, u;

    block_file->file = (char *) NULL;
    block_file->fp = (FILE *) NULL;

    STRING_MALLOC_COPY(block_file->file, file);

    COPY_VECTOR(block_file->points, points, ndim);
    COPY_VECTOR(block_file->block_size, block_size, ndim);

    BLOCKS(block_file->blocks, points, block_file->block_size, ndim);

    CUMULATIVE(block_file->cum_points, points, block_file->npoints, ndim);
    CUMULATIVE(block_file->cum_blocks, block_file->blocks, block_file->nblocks, ndim);

    VECTOR_PRODUCT(u, block_size, ndim);

    t = 1;
    for (i = 0; i < ndim; i++)
    {
	if (padded)
	    s = block_size[i] * block_file->blocks[i];
	else
	    s = points[i];

	block_file->cum_block_points[i] = u * t;

	t *= s;
	u /= block_size[i];
    }

    return CCPN_OK;
}
	
Block_file new_block_file(char *file, int ndim,
	int *points, int *block_size, Bool *dim_wrapped,
	Mem_cache mem_cache, int bytes_per_point, Bool big_endian,
	Bool padded, int header, Bool integer, Bool writeable, int block_header)
{
    Block_file block_file = (Block_file) NULL;

    if ((ndim > MAX_NDIM) || (bytes_per_point > BYTES_PER_WORD))
	return NULL;

    if (writeable && bytes_per_point != BYTES_PER_WORD)
	return NULL;

    MALLOC_NEW(block_file, struct Block_file, 1);
    block_file->block_file_kind = BLOCK_FILE;

    if (init_block_file(block_file, file, ndim,
				points, block_size, padded) == CCPN_ERROR)
    {
	delete_block_file(block_file);
	return NULL;
    }

    block_file->ndim = ndim;
    block_file->mem_cache = mem_cache;
    block_file->bytes_per_point = bytes_per_point;
    block_file->big_endian = big_endian;
    block_file->padded = padded;
    block_file->header = header;
    block_file->integer = integer;
    block_file->writeable = writeable;
    block_file->block_header = block_header;
    COPY_VECTOR(block_file->dim_wrapped, dim_wrapped, ndim);

    block_file->block_table = new_hash_table(equal_block, hash_block);
    block_file->dirty_blocks = new_hash_table(equal_block, hash_block);

    if (!block_file->block_table)
    {
	delete_block_file(block_file);
	return NULL;
    }

    block_file->fp = (FILE *) NULL;

    return block_file;
}

Block_file new_shape_block_file(char *file, int ndim, int *points,
	int *block_size, Bool *dim_wrapped, Shape_file shape_file)
{
/* TBD: should be FALSE but something not working then
    Bool padded = CCPN_FALSE;
*/
    Bool padded = CCPN_TRUE;
    Shape_block_file shape_block_file;
    Block_file block_file;

    if (ndim > MAX_NDIM)
	return NULL;

    MALLOC_NEW(shape_block_file, struct Shape_block_file, 1);
    block_file = (Block_file) shape_block_file;
    block_file->block_file_kind = SHAPE_FILE;
    shape_block_file->shape_file = NULL;

    if (init_block_file(block_file, file, ndim,
			points, block_size, padded) == CCPN_ERROR)
    {
	delete_block_file(block_file);
	return NULL;
    }

    shape_block_file->shape_file = shape_file;
    block_file->ndim = ndim;
    block_file->padded = padded;
    COPY_VECTOR(block_file->dim_wrapped, dim_wrapped, ndim);

    return block_file;
}

static void remove_block(Hash_key key, Hash_data data, Generic_ptr user_data)
{
    Mem_cache mem_cache = (Mem_cache) user_data;
    Block_data block_data = (Block_data) data;

/*  below automatically calls delete_data_block  */
    remove_mem_cache(mem_cache, (Generic_ptr) block_data);
}

void delete_block_file(Block_file block_file)
{
/*
    printf("delete_block_file\n");
*/

    if (block_file && block_file->block_file_kind == BLOCK_FILE)
    {
	close_block_file(block_file);
	clear_hash_table(block_file->block_table, remove_block, (Generic_ptr) block_file->mem_cache);
	delete_hash_table(block_file->block_table);
	delete_hash_table(block_file->dirty_blocks);
    }
    else if (block_file)
    {
	delete_shape_file(((Shape_block_file) block_file)->shape_file);
    }

    FREE(block_file->file, char);
    FREE(block_file, struct Block_file);
}

CcpnStatus open_block_file(Block_file block_file)
{
    Integer m, n;
#define  MAX_ZEROES  1024
    char zeroes[MAX_ZEROES];
    Bool zero_file;

    if (block_file->block_file_kind == SHAPE_FILE)
	return CCPN_OK;

    close_block_file(block_file);

    if (block_file->writeable)
    {
/* always zero file for now
	if (file_exists(block_file->file))
	    zero_file = CCPN_FALSE;
	else
*/
	    zero_file = CCPN_TRUE;

	if (OPEN_FOR_BINARY_WRITE_UPDATE(block_file->fp, block_file->file))
	    return CCPN_ERROR;

	if (zero_file)
	{
	    n = block_file->bytes_per_point * block_file->npoints;
	    n += block_file->header;
	    ZERO_VECTOR(zeroes, MAX_ZEROES);
	    for (m = 0; n > 0; n -= m)
	    {
		m = MIN(n, MAX_ZEROES);
		FWRITE(zeroes, 1, m, block_file->fp);
	    }
	    rewind(block_file->fp);
	}
    }
    else
    {
	if (OPEN_FOR_BINARY_READING(block_file->fp, block_file->file))
	    return CCPN_ERROR;
    }

    return CCPN_OK;
}

void close_block_file(Block_file block_file)
{
    if (block_file->block_file_kind == SHAPE_FILE)
	return;

    FCLOSE(block_file->fp);
}

CcpnStatus have_peak_block_file(Block_file block_file,
	int *first, int *last, Bool have_low, Bool have_high,
	float low, float high, Bool *have_peak, CcpnString error_msg)
{
    int i, npoints, ndim= block_file->ndim;
    int point[MAX_NDIM];
    Integer cum_points[MAX_NDIM];
    float v;

    *have_peak = CCPN_FALSE;

    if (!have_low && !have_high)
	return CCPN_OK;

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {
	cum_points[i] = npoints;
	npoints *= last[i] - first[i];
    }
 
    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, first, ndim);

	CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));
 
	if ((have_high && (v >= high)) || (have_low && (v <= low)))
	{
	    *have_peak = CCPN_TRUE;
	    break;
	}
    }

    return CCPN_OK;
}

CcpnStatus max_value_block_file(Block_file block_file,
	int *first, int *last, float *max_value, CcpnString error_msg)
{
    int i, npoints, ndim= block_file->ndim;
    int point[MAX_NDIM];
    Integer cum_points[MAX_NDIM];
    float v;

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {
        if ((first[i] < 0) || (first[i] > block_file->points[i]))
        {
            sprintf(error_msg, "first[%d] = %d, must be between 0 and %d",
        			i, first[i], block_file->points[i]);
            return CCPN_ERROR;
        }

        if ((last[i] < 0) || (last[i] > block_file->points[i]))
        {
            sprintf(error_msg, "last[%d] = %d, must be between 0 and %d",
        			i, last[i], block_file->points[i]);
            return CCPN_ERROR;
        }

        if (first[i] >= last[i])
        {
            sprintf(error_msg, "first[%d] = %d >= last[%d] = %d",
        			i, first[i], i, last[i]);
            return CCPN_ERROR;
        }

        cum_points[i] = npoints;
        npoints *= last[i] - first[i];

    }

 
    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, first, ndim);

	CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));

        *max_value = MAX(*max_value, v);
    }

    return CCPN_OK;
}

CcpnStatus min_value_block_file(Block_file block_file,
	int *first, int *last, float *min_value, CcpnString error_msg)
{
    int i, npoints, ndim= block_file->ndim;
    int point[MAX_NDIM];
    Integer cum_points[MAX_NDIM];
    float v;

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {

        if ((first[i] < 0) || (first[i] > block_file->points[i]))
        {
            sprintf(error_msg, "first[%d] = %d, must be between 0 and %d",
        			i, first[i], block_file->points[i]);
            return CCPN_ERROR;
        }

        if ((last[i] < 0) || (last[i] > block_file->points[i]))
        {
            sprintf(error_msg, "last[%d] = %d, must be between 0 and %d",
        			i, last[i], block_file->points[i]);
            return CCPN_ERROR;
        }

        if (first[i] >= last[i])
        {
            sprintf(error_msg, "first[%d] = %d >= last[%d] = %d",
        			i, first[i], i, last[i]);
            return CCPN_ERROR;
        }

        cum_points[i] = npoints;
        npoints *= last[i] - first[i];
    }
 
    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, first, ndim);

	CHECK_STATUS(get_point_block_file(block_file, &v, point, error_msg));

        *min_value = MIN(*min_value, v);
    }

    return CCPN_OK;
}

/* this locks returned block_data */
Block_data get_block_data_number(Block_file block_file, Integer block)
{
    int block_array[MAX_NDIM];

    ARRAY_OF_INDEX(block_array, block, block_file->cum_blocks, block_file->ndim);

    return get_block_data_array(block_file, block_array);
}

static Block_data get_shape_block_data(Shape_block_file shape_block_file,
	int *block_array, int ncomponents, int *components)
{
    Block_file block_file = (Block_file) shape_block_file;
    int i, j, npts, ndim = block_file->ndim;
    int first[MAX_NDIM], last[MAX_NDIM], array[MAX_NDIM];
    int *block_size = block_file->block_size;
    int *points = block_file->points;
    Integer *cum_block_size;
    float value;
    Int_array block;
    Block_data block_data;
    Bool flag;

    for (i = 0; i < ndim; i++)
    {
        first[i] = block_array[i] * block_size[i];
        last[i] = first[i] + block_size[i];
        last[i] = MIN(last[i], points[i]);
        if (first[i] >= last[i])
            return (Block_data) NULL;
    }

    block = new_int_array(ndim, block_array);

    if (!block)
	return (Block_data) NULL;

    block_data = new_block_data(block_file, block);

    if (!block_data)
        return (Block_data) NULL;

    cum_block_size = block_data->cum_block_size;
    VECTOR_PRODUCT(npts, block_size, ndim);

    for (i = 0; i < npts; i++)
    {
        ARRAY_OF_INDEX(array, i, cum_block_size, ndim);
        ADD_VECTORS(array, array, first, ndim);
/* TBD: have flag here only because now assuming padded = TRUE */
        flag = CCPN_FALSE;
        for (j = 0; j < ndim; j++)
	{
	    if (array[j] >= points[j])
	    {
		flag = CCPN_TRUE;
		break;
	    }
	}
	if (flag)
	    value = 0;
        else
            value = get_point_shape_file(shape_block_file->shape_file, array, ncomponents, components);
        block_data->data[i] = value;
    }

    return block_data;
}

/* this locks returned block_data */
Block_data get_block_data_array(Block_file block_file, int *block_array)
{
    int i;
    Hash_data data;
    Block_data block_data;
    Int_array block;

    if (block_file->block_file_kind == SHAPE_FILE)
	return get_shape_block_data((Shape_block_file) block_file, block_array, -1, NULL);

    for (i = 0; i < block_file->ndim; i++)
    {
	if ((block_array[i] < 0) || (block_array[i] >= block_file->blocks[i]))
	    return (Block_data) NULL;
    }

    block = new_int_array(block_file->ndim, block_array);

    if (!block)
	return (Block_data) NULL;

    if (find_hash_table(block_file->block_table, (Hash_key) block, &data))
    {
	block_data = (Block_data) data;
	if (lock_block_data(block_file, block_data) == CCPN_ERROR)
	    return (Block_data) NULL;
    }

    /* do it again in case deleted before lock obtained */
    if (data &&
	find_hash_table(block_file->block_table, (Hash_key) block, &data))
    {
	block_data = (Block_data) data;
	delete_int_array(block);
    }
    else
    {
	block_data = new_block_data(block_file, block);

	if (!block_data)
	    return (Block_data) NULL;

	if (read_block_data(block_file, block_data) == CCPN_ERROR)
	{
	    delete_block_data(block_data);
	    return (Block_data) NULL;
	}

        if (insert_hash_table(block_file->block_table,
		(Hash_key) block_data->block, (Hash_data) block_data) == CCPN_ERROR)
	{
	    delete_block_data(block_data);
	    return (Block_data) NULL;
	}

	if (add_mem_cache(block_file->mem_cache,
		(Generic_ptr) block_data, sizeof(Float)*block_data->size,
		delete_data_block, block_file) == CCPN_ERROR)
	{
	    delete_data_block((Generic_ptr) block_data, (Generic_ptr) block_file);
	    return (Block_data) NULL;
	}
    }

    return block_data;
}

Block_data get_components_block_data(Block_file block_file, int *block_array,
	int ncomponents, int *components)
{
    if (block_file->block_file_kind == BLOCK_FILE)
	return get_block_data_array(block_file, block_array);
    else
	return get_shape_block_data((Shape_block_file) block_file,
				block_array, ncomponents, components);
}

CcpnStatus lock_block_data(Block_file block_file, Block_data block_data)
{
    return lock_mem_cache(block_file->mem_cache, (Generic_ptr) block_data);
}

CcpnStatus unlock_block_data(Block_file block_file, Block_data block_data)
{
    return unlock_mem_cache(block_file->mem_cache, (Generic_ptr) block_data);
}

CcpnStatus get_point_block_file(Block_file block_file, float *v, int *point,
				CcpnString error_msg)
{
    int i, p, ndim = block_file->ndim, block[MAX_NDIM], remainder[MAX_NDIM];
    Block_data block_data;
    Line msg;

    for (i = 0; i < ndim; i++)
    {
	if ((point[i] < 0) || (point[i] >= block_file->points[i]))
        {
	    sprintf(error_msg, "getting point %d in dimension %d: out of range (1, %d)",
			point[i]+1, i+1, block_file->points[i]);
	    return CCPN_ERROR;
        }
    }

    for (i = 0; i < ndim; i++)
	block[i] = point[i] /  block_file->block_size[i];

    block_data = get_block_data_array(block_file, block);

    if (!block_data)
    {
        sprintf(error_msg, "could not get point:");
        for (i = 0; i < ndim; i++)
        {
            sprintf(msg, " %d", point[i]+1);
            strcat(error_msg, msg);
        }

	return CCPN_ERROR;
    }

    for (i = 0; i < ndim; i++)
	remainder[i] = point[i] % block_file->block_size[i];

    INDEX_OF_ARRAY(p, remainder, block_data->cum_block_size, ndim);

    *v = block_data->data[p];

    if (block_file->block_file_kind == BLOCK_FILE)
	(void) unlock_block_data(block_file, block_data);

    return CCPN_OK;
}

static CcpnStatus get_shape_point_block_file(Shape_block_file shape_block_file,
	float *v, int *point, int ncomponents, int *components, CcpnString error_msg)
{
    Block_file block_file = (Block_file) shape_block_file;
    int i, p, ndim = block_file->ndim, block[MAX_NDIM], remainder[MAX_NDIM];
    Block_data block_data;
    Line msg;

    for (i = 0; i < ndim; i++)
    {
	if ((point[i] < 0) || (point[i] >= block_file->points[i]))
        {
	    sprintf(error_msg, "getting point %d in dimension %d: out of range (1, %d)",
			point[i]+1, i+1, block_file->points[i]);
	    return CCPN_ERROR;
        }
    }

    for (i = 0; i < ndim; i++)
	block[i] = point[i] /  block_file->block_size[i];

    block_data = get_shape_block_data(shape_block_file, block, ncomponents, components);

    if (!block_data)
    {
        sprintf(error_msg, "could not get point:");
        for (i = 0; i < ndim; i++)
        {
            sprintf(msg, " %d", point[i]+1);
            strcat(error_msg, msg);
        }

	return CCPN_ERROR;
    }

    for (i = 0; i < ndim; i++)
	remainder[i] = point[i] % block_file->block_size[i];

    INDEX_OF_ARRAY(p, remainder, block_data->cum_block_size, ndim);

    *v = block_data->data[p];

    return CCPN_OK;
}

CcpnStatus get_components_point_block_file(Block_file block_file, float *v, int *point,
			int ncomponents, int *components, CcpnString error_msg)
{
    if (block_file->block_file_kind == BLOCK_FILE)
	return get_point_block_file(block_file, v, point, error_msg);
    else
	return get_shape_point_block_file((Shape_block_file) block_file,
			v, point, ncomponents, components, error_msg);
}

CcpnStatus get_position_block_file(Block_file block_file, float *v,
					float *position, CcpnString error_msg)
{
    int i, j, npoints, ndim= block_file->ndim;
    int point[MAX_NDIM], first[MAX_NDIM], pnt[MAX_NDIM];
    Integer cum_points[MAX_NDIM];
    float value, z, w, diff[MAX_NDIM];

    for (i = 0; i < ndim; i++)
    {
	if ((position[i] < 0) || (position[i] >= block_file->points[i]))
	{
	    sprintf(error_msg, "position[%d] = %4.3f, must be in [0, %d)",
				i, position[i], block_file->points[i]);
	    return CCPN_ERROR;
	}
    }

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {
	cum_points[i] = npoints;
	npoints *= 2;

	first[i] = position[i];
    }
 
/*
printf("position = %f %f %f\n", position[0], position[1], position[2]);
*/
    value = 0;
    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, first, ndim);

	for (j = 0; j < ndim; j++)
	    pnt[j] = point[j] % block_file->points[j];

	CHECK_STATUS(get_point_block_file(block_file, &z, pnt, error_msg));

	SUBTRACT_VECTORS(diff, point, position, ndim);
	w = 1.0;
        for (j = 0; j < ndim; j++)
	    w *= 1.0 - ABS(diff[j]);

	value = value + w*z;
/*
printf("i = %d, point = %d %d %d\n", i, point[0], point[1], point[2]);
printf("i = %d, diff = %f %f %f\n", i, diff[0], diff[1], diff[2]);
printf("i = %d, w = %f, z = %f, value = %f\n", i, w, z, value);
*/
    }

    *v = value;

    return CCPN_OK;
}

static CcpnStatus half_max_position(Block_file block_file, Bool find_maximum,
		float v, int *point, int dim, int dirn,
		float *half_max, CcpnString error_msg)
{
    int i, i_start, i_end, i_step;
    float v_half = 0.5 * v, v_prev = v, v_this;
    int q[MAX_NDIM];

/*
printf("half_max_position1: v = %3.2f, dim = %d, dirn = %d, point = %d %d %d\n",
v, dim, dirn, point[0], point[1], point[2]);
*/
    if (dirn == 1)
    {
	i_start = point[dim] + 1;
	if (block_file->dim_wrapped[dim])
	    i_end = i_start + block_file->points[dim] - 1;
	else
	    i_end = block_file->points[dim];
	i_step = 1;
    }
    else
    {
	i_start = point[dim] - 1;
	if (block_file->dim_wrapped[dim])
	    i_end = i_start - block_file->points[dim] + 1;
	else
	    i_end = -1;
	i_step = -1;
    }
/*
printf("half_max_position2: i_start = %d, i_end = %d, i_step = %d\n",
i_start, i_end, i_step);
*/

    COPY_VECTOR(q, point, block_file->ndim);

    for (i = i_start; i != i_end; i += i_step)
    {
	q[dim] = (i + block_file->points[dim]) % block_file->points[dim];
/*
printf("half_max_position3: i = %d, q = %d %d %d\n", i, q[0], q[1], q[2]);
*/
	CHECK_STATUS(get_point_block_file(block_file, &v_this, q, error_msg));
/*
printf("half_max_position4: i = %d, v_this = %3.2f, v_prev = %3.2f\n",
i, v_this, v_prev);
*/

	if (find_maximum)
	{
	    if (v_this < v_half)
	    {
		*half_max = i - i_step*(v_half-v_this)/(v_prev-v_this);
		return CCPN_OK;
	    }
	}
	else
	{
	    if (v_this > v_half)
	    {
		*half_max = i - i_step*(v_half-v_this)/(v_prev-v_this);
		return CCPN_OK;
	    }
	}

	v_prev = v_this;
    }

    if (dirn == 1)
	*half_max = block_file->points[i] - 1.0;
    else
	*half_max = 1.0;

    return CCPN_OK;
}

CcpnStatus linewidth_block_file(Block_file block_file,
			Bool have_maximum, float v, int *point, int dim,
			float *linewidth, CcpnString error_msg)
{
    float a, b;

    if (half_max_position(block_file, have_maximum, v, point, dim, 1,
				&a, error_msg) == CCPN_ERROR)
	return CCPN_ERROR;

    if (half_max_position(block_file, have_maximum, v, point, dim, -1,
				&b, error_msg) == CCPN_ERROR)
	return CCPN_ERROR;

    *linewidth = a - b;

    return CCPN_OK;
}

CcpnStatus get_box_block_file(Block_file block_file, int *p_npoints, float **p_v,
                        int *box_min, int *box_max, CcpnString error_msg)
{
    int i, npoints, ndim= block_file->ndim;
    int point[MAX_NDIM];
    Integer cum_points[MAX_NDIM];
    float z, *v;

    for (i = 0; i < ndim; i++)
    {
	if ((box_min[i] < 0) || (box_min[i] > block_file->points[i]))
	{
	    sprintf(error_msg, "box_min[%d] = %d, must be between 0 and %d",
				i, box_min[i], block_file->points[i]);
	    return CCPN_ERROR;
	}

	if ((box_max[i] < 0) || (box_max[i] > block_file->points[i]))
	{
	    sprintf(error_msg, "box_max[%d] = %d, must be between 0 and %d",
				i, box_max[i], block_file->points[i]);
	    return CCPN_ERROR;
	}

	if (box_min[i] > box_max[i])
	{
	    sprintf(error_msg, "box_min[%d] = %d > box_max[%d] = %d",
				i, box_min[i], i, box_max[i]);
	    return CCPN_ERROR;
	}
    }

    npoints = 1;
    for (i = 0; i < ndim; i++)
    {
	cum_points[i] = npoints;
	npoints *= box_max[i] - box_min[i];
    }
 
    *p_v = NULL;
    *p_npoints = npoints;
    if (npoints == 0)
      return CCPN_OK;

    sprintf(error_msg, "allocating box memory");
    MALLOC(v, float, npoints);
    *p_v = v;

    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, box_min, ndim);

	if (get_point_block_file(block_file, &z, point, error_msg) == CCPN_ERROR)
	{
	    FREE(v, float);
	    return CCPN_ERROR;
	}

        v[i] = z;
    }

    return CCPN_OK;
}

CcpnStatus set_point_block_file(Block_file block_file, float v, int *point, CcpnString error_msg)
{
    int i, p, ndim = block_file->ndim, block[MAX_NDIM], remainder[MAX_NDIM];
    Block_data block_data;
    Line msg;

    sprintf(error_msg, "setting point: ");
    for (i = 0; i < ndim; i++)
    {
        sprintf(msg, " %d", point[i]+1);
        strcat(error_msg, msg);
    }

    for (i = 0; i < ndim; i++)
	block[i] = point[i] /  block_file->block_size[i];

    block_data = get_block_data_array(block_file, block);

    if (!block_data)
    {
	strcat(error_msg, ": getting block_data");
	return CCPN_ERROR;
    }

    for (i = 0; i < ndim; i++)
	remainder[i] = point[i] % block_file->block_size[i];

    INDEX_OF_ARRAY(p, remainder, block_data->cum_block_size, ndim);

    if (contains_hash_table(block_file->dirty_blocks, (Hash_key) block_data->block))
    {
	/* decrement count because only want one link to it */
	(void) unlock_block_data(block_file, block_data);
    }
    else
    {
	strcat(error_msg, ": inserting into hash_table");
	CHECK_STATUS(insert_hash_table(block_file->dirty_blocks, (Hash_key) block_data->block, (Hash_data) block_data));
    }

    block_data->data[p] = v;

    return CCPN_OK;
}

CcpnStatus set_box_block_file(Block_file block_file, int npoints, float *v,
				int *box_min, int *box_max, CcpnString error_msg)
{
    int i, j, npts, ndim= block_file->ndim;
    int point[MAX_NDIM];
    Integer cum_points[MAX_NDIM];

    for (i = 0; i < ndim; i++)
    {
	if ((box_min[i] < 0) || (box_min[i] > block_file->points[i]))
	{
	    sprintf(error_msg, "box_min[%d] = %d, must be between 0 and %d",
				i, box_min[i], block_file->points[i]);
	    return CCPN_ERROR;
	}

	if ((box_max[i] < 0) || (box_max[i] > block_file->points[i]))
	{
	    sprintf(error_msg, "box_max[%d] = %d, must be between 0 and %d",
				i, box_max[i], block_file->points[i]);
	    return CCPN_ERROR;
	}

	if (box_min[i] > box_max[i])
	{
	    sprintf(error_msg, "box_min[%d] = %d > box_max[%d] = %d",
				i, box_min[i], i, box_max[i]);
	    return CCPN_ERROR;
	}
    }

    npts = 1;
    for (i = 0; i < ndim; i++)
    {
	cum_points[i] = npts;
	npts *= box_max[i] - box_min[i];
    }
 
    if (npts != npoints)
    {
	sprintf(error_msg, "npoints = %d, inconsistent with box size = %d",
							npoints, npts);
	return CCPN_ERROR;
    }

    for (i = 0; i < npoints; i++)
    {
	ARRAY_OF_INDEX(point, i, cum_points, ndim);
	ADD_VECTORS(point, point, box_min, ndim);

	CHECK_STATUS(set_point_block_file(block_file, v[i], point, error_msg));
    }

    return CCPN_OK;
}

static CcpnStatus save_block(Hash_key key, Hash_data data, Generic_ptr user_data, CcpnString error_msg)
{
    int block;
    CcpnStatus status;
    Block_data block_data = (Block_data) data;
    Block_file block_file = (Block_file) user_data;
 
    status = write_block_data(block_file, block_data);

    if (status == CCPN_ERROR)
    {
	INDEX_OF_ARRAY(block, block_data->block->values, block_file->cum_blocks, block_file->ndim);
	sprintf(error_msg, "writing block %d", block);
    }

    return status;
}

static void clear_block(Hash_key key, Hash_data data, Generic_ptr user_data)
{
    Block_data block_data = (Block_data) data;
    Block_file block_file = (Block_file) user_data;

    (void) unlock_block_data(block_file, block_data);
}

CcpnStatus save_block_file(Block_file block_file, CcpnString error_msg)
{
    if (!block_file->writeable)
	RETURN_ERROR_MSG("block_file not writeable");

    CHECK_STATUS(iterate_hash_table(block_file->dirty_blocks, save_block, (Generic_ptr) block_file, error_msg));
    clear_hash_table(block_file->dirty_blocks, clear_block, (Generic_ptr) block_file);

    return CCPN_OK;
}

static CcpnStatus check_func(Hash_key key, Hash_data data, Generic_ptr user_data, CcpnString error_msg)
{
    int i;
    Block_data block_data = (Block_data) data;

    printf("block_file check_func: %p:", block_data);
    for (i = 0; i < block_data->ndim; i++)
	printf(" %d", block_data->block->values[i]);
    printf("\n");

    return CCPN_OK;
}

CcpnStatus check_block_file(Block_file block_file, CcpnString error_msg)
{
    printf("check_block_file\n");
    CHECK_STATUS(iterate_hash_table(block_file->block_table, check_func, NULL, error_msg));
    CHECK_STATUS(check_mem_cache(block_file->mem_cache, error_msg));

    return CCPN_OK;
}

static CcpnStatus get_slice(Block_file block_file, int dim, int *slice,
				int *p_npoints, float **p_w, CcpnString error_msg)
{
    int i, j, box_min[MAX_NDIM], box_max[MAX_NDIM];

    for (i = j = 0; i < block_file->ndim; i++)
    {
	if (i == dim)
	{
	    box_min[i] = 0;
	    box_max[i] = block_file->points[i];
	}
	else
	{
	    box_min[i] = slice[j++] - 1;
	    box_max[i] = slice[j++];
	}
    }

    return get_box_block_file(block_file, p_npoints, p_w, box_min, box_max, error_msg);
}

CcpnStatus compare_slices_block_files(Block_file block_file1, int dim1, int *slice1,
 		Block_file block_file2, int dim2, int *slice2, float *v,
		float scale, float offset, CcpnString error_msg)
{
    int i, j, npoints1, npoints2;
    float a, b, c, d, x1, x2, *w1, *w2;

    if ((dim1 < 0) || (dim1 >= block_file1->ndim))
    {
        sprintf(error_msg, "dim1 = %d, should be >= 0 and < %d", dim1, block_file1->ndim);
        return CCPN_ERROR;
    }

    if ((dim2 < 0) || (dim2 >= block_file2->ndim))
    {
        sprintf(error_msg, "dim2 = %d, should be >= 0 and < %d", dim2, block_file2->ndim);
        return CCPN_ERROR;
    }

    if (get_slice(block_file1, dim1, slice1, &npoints1, &w1, error_msg) == CCPN_ERROR)
	return CCPN_ERROR;

    if (get_slice(block_file2, dim2, slice2, &npoints2, &w2, error_msg) == CCPN_ERROR)
    {
	FREE(w1, float);
	return CCPN_ERROR;
    }

    a = b = c = 0;
    for (i = 0; i < npoints1; i++)
    {
	d = scale*i + offset;
	if ((d >= 0) && (d <= npoints2-1))
	{
	    j = (int) d;
	    d -= j;
	    x1 = w1[i];
	    x2 = w2[j]*(1-d) + w2[j+1]*d;
	    a += x1 * x1;
	    b += x2 * x2;
	    c += x1 * x2;
	}
    }

    FREE(w1, float);
    FREE(w2, float);

    if (a > 0 && b > 0)
    {
    	*v = c / (float) sqrt((double) (a*b));
    }
    else
    {
	*v = 0;
    }

    return CCPN_OK;
}

Block_plane *alloc_block_planes(int nplanes)
{
    Block_plane *block_planes;

    MALLOC_NEW(block_planes, Block_plane, nplanes);

    return block_planes;
}

void free_block_planes(Block_plane *block_planes)
{
    FREE(block_planes, Block_plane);
}

CcpnStatus plane_product_block_files(Block_file block_file,
	int nplanes, Block_plane *block_planes, float scale, CcpnString error_msg)
{
    int nx, ny, i, j, k, m, dim, point[MAX_NDIM];
    float prod, maxprod = 0, v, pst, position[MAX_NDIM], *planeProd;
    Bool bad;
    Block_plane *bp;
    Block_file bf;

    if (!block_file->writeable)
	RETURN_ERROR_MSG("output block_file not writeable");

    if (block_file->ndim != 2)
	RETURN_ERROR_MSG("output block_file must be 2D");

    nx = block_file->points[0];
    ny = block_file->points[1];

    if (nx < 2)
	RETURN_ERROR_MSG("output block_file must have at least 2 points in xdim");

    if (ny < 2)
	RETURN_ERROR_MSG("output block_file must have at least 2 points in ydim");

    sprintf(error_msg, "allocating planeProd");
    MALLOC(planeProd, float, nx*ny);

    for (j = 0; j < ny; j++)
    {
	for (i = 0; i < nx; i++)
	{
	    prod = 1.0;
	    bad = CCPN_FALSE;
	    for (k = 0; k < nplanes; k++)
	    {
		bp = block_planes + k;
		bf = bp->block_file;
		m = 0;
		for (dim = 0; dim < bf->ndim; dim++)
		{
		    if (dim == bp->xdim)
			pst = bp->xmin + i*(bp->xmax - bp->xmin)/(nx-1);
		    else if (dim == bp->ydim)
			pst = bp->ymin + j*(bp->ymax - bp->ymin)/(ny-1);
		    else
			pst = bp->plane[m++];

		    if ((pst < 0) || (pst > bf->points[dim]-1))
		    {
			bad = CCPN_TRUE;
			prod = 0.0;
			break;
		    }

		    position[dim] = pst;
		}

		if (bad)
		    break;

		if (get_position_block_file(bf, &v, position, error_msg) == CCPN_ERROR)
		{
		    FREE(planeProd, float);
		    return CCPN_ERROR;
		}

		prod *= v;
	    }

	    maxprod = MAX(maxprod, ABS(prod));
	    planeProd[nx*j+i] = prod;
	}
    }

    if (maxprod > 0)  // if not irrelevant because all data is 0
	scale /= maxprod;
    else
	scale = 1;


    for (j = 0; j < ny; j++)
    {
	point[1] = j;
	for (i = 0; i < nx; i++)
	{
	    point[0] = i;
	    prod = scale * planeProd[nx*j+i];
	    if (set_point_block_file(block_file, prod, point, error_msg) == CCPN_ERROR)
	    {
		FREE(planeProd, float);
		return CCPN_ERROR;
	    }
	}
    }

    FREE(planeProd, float);

    CHECK_STATUS(save_block_file(block_file, error_msg));

    return CCPN_OK;
}


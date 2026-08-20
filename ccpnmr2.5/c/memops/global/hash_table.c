
/*
======================COPYRIGHT/LICENSE START==========================

hash_table.c: Part of the CcpNmr Analysis program

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
#include "hash_table.h"

#define TABLE_INITIAL_SIZE	931
#define GROW_FACTOR		2
#define RESIZE_FACTOR		3
#define SHRINK_FACTOR		64

#define  HASH_MASK		((unsigned long) (1<<30) - 1)

typedef struct Hash_entry
{
    Bool used;
    Hash_key key;
    Hash_data data;
    Hash_value hash;
}   Hash_entry;

/* not great hash functions but use for now... */

static Hash_value combine_hash(Hash_value hash1, Hash_value hash2)
{
/*  HASH_MASK used so that arithmetic guaranteed not to overflow  */
/*  lose two bits (could lose only one but safer this way, because of pos/neg)  */
    return ((hash1&HASH_MASK) + 123456789) ^ hash2;
}

static Hash_value create_hash(int value)
{
    return value ^ (value << 3) ^ (value << 7);
}

Bool equal_pointers(Hash_key key1, Hash_key key2)
{
    return key1 == key2;
}

Hash_value hash_pointers(Hash_key key)
{
    return create_hash((int) key);
}

static Hash_entry *next_entry(Hash_table hash_table, Hash_entry *entry)
{
    if ((entry + 1) < (hash_table->entries + hash_table->nentries))
	return (entry + 1);
    else
	return hash_table->entries;
}

static Hash_entry *table_entry(Hash_table hash_table, Hash_key key,
							Hash_value hash)
{
    Hash_entry *entry;

    if (hash_table->nentries == 0)
	return NULL;

    for (entry = hash_table->entries + (hash % hash_table->nentries);
			entry->used; entry = next_entry(hash_table, entry))
    {
	if ((entry->hash == hash) && (hash_table->equal_func)(entry->key, key))
	    return entry;
    }

    return entry;
}

static void insert_table(Hash_table hash_table, Hash_key key, Hash_data data, Hash_value hash)
{
    Hash_entry *entry = table_entry(hash_table, key, hash);

/*
    printf("insert_table1: nentries = %d, nused = %d\n",
			hash_table->nentries, hash_table->nused);
*/
    if (!entry->used)
    {
	entry->used = CCPN_TRUE;
	hash_table->nused++;
    }

    entry->key = key;
    entry->data = data;
    entry->hash = hash;

/*
    printf("insert_table2: nentries = %d, nused = %d\n",
			hash_table->nentries, hash_table->nused);
*/
}

static CcpnStatus resize_table(Hash_table hash_table)
{
    int i, nentries;
    Hash_entry *old_entries = hash_table->entries, *new_entries;

/*
    printf("resize_table1: nentries = %d, nused = %d\n",
			hash_table->nentries, hash_table->nused);
*/
    if (hash_table->nentries == 0)
	nentries = TABLE_INITIAL_SIZE;
    else
	nentries = TABLE_INITIAL_SIZE + RESIZE_FACTOR * hash_table->nused;

    MALLOC(new_entries, Hash_entry, nentries);

    for (i = 0; i < nentries; i++)
	new_entries[i].used = CCPN_FALSE;

    hash_table->nused = 0;
    hash_table->entries = new_entries;
    SWAP(hash_table->nentries, nentries, int);

    for (i = 0; i < nentries; i++)
    {
	if (old_entries[i].used)
	    insert_table(hash_table, old_entries[i].key,
				old_entries[i].data, old_entries[i].hash);
    }

    FREE(old_entries, Hash_entry);

/*
    printf("resize_table2: nentries = %d, nused = %d\n",
			hash_table->nentries, hash_table->nused);

    for (i = 0; i < hash_table->nentries; i++)
    {
        if (hash_table->entries[i].used)
            printf("resize_table3: i = %d, used = %d, key = %d, data = %d\n",
                        i, hash_table->entries[i].used,
			(int) (hash_table->entries[i].key),
			(int) (hash_table->entries[i].data));
    }
*/

    return CCPN_OK;
}

Hash_table new_hash_table(Equal_func equal_func, Hash_func hash_func)
{
    Hash_table hash_table = (Hash_table) NULL;

    MALLOC_NEW(hash_table, struct Hash_table, 1);

    if (equal_func)
	hash_table->equal_func = equal_func;
    else
	hash_table->equal_func = equal_pointers;

    if (hash_func)
	hash_table->hash_func = hash_func;
    else
	hash_table->hash_func = hash_pointers;

    hash_table->nentries = 0;
    hash_table->nused = 0;
    hash_table->entries = (Hash_entry *) NULL;

    return hash_table;
}

void delete_hash_table(Hash_table hash_table)
{
    if (hash_table)
    {
	FREE(hash_table->entries, Hash_entry);
	FREE(hash_table, struct Hash_table);
    }
}

CcpnStatus insert_hash_table(Hash_table hash_table, Hash_key key, Hash_data data)
{
/*
    int i;
*/
    Hash_value hash;

    if (hash_table->nused >= hash_table->nentries / GROW_FACTOR)
	CHECK_STATUS(resize_table(hash_table));

    hash = (hash_table->hash_func)(key);

    insert_table(hash_table, key, data, hash);

/*
    for (i = 0; i < hash_table->nentries; i++)
    {
        if (hash_table->entries[i].used)
            printf("insert_hash_table: i = %d, key = %d, data = %d\n",
                        i, (int) (hash_table->entries[i].key),
                        (int) (hash_table->entries[i].data));
    }
*/

    return CCPN_OK;
}

Hash_data remove_hash_table(Hash_table hash_table, Hash_key key)
{
/*
    int i;
*/
    Hash_value hash = (hash_table->hash_func)(key);
    Hash_entry *entry = table_entry(hash_table, key, hash), *ne, *pe;
    Hash_data data = (Hash_data) NULL;

/*
    printf("remove_hash_table0: hash = %ld, entry->key = %d, entry->data = %d\n",
                        hash, (int) (entry->key), (int) (entry->data));

    for (i = 0; i < hash_table->nentries; i++)
    {
        if (hash_table->entries[i].used)
            printf("remove_hash_table1: i = %d, key = %d, data = %d\n",
                        i, (int) (hash_table->entries[i].key),
                        (int) (hash_table->entries[i].data));
    }
*/

    if (entry && entry->used)
    {
	data = entry->data;
	entry->used = CCPN_FALSE;

	for (ne = next_entry(hash_table, entry); ne->used;
					ne = next_entry(hash_table, ne))
	{
	    pe = hash_table->entries + (ne->hash % hash_table->nentries);
	    if ((entry < ne) ? ((pe <= entry) || (pe > ne))
	    		     : ((pe <= entry) && (pe > ne)))
	    {
		*entry = *ne;
		entry->used = CCPN_TRUE;
		entry->key = ne->key;
		entry->data = ne->data;
		entry->hash = ne->hash;
		ne->used = CCPN_FALSE;
		entry = ne;
	    }
	}

	hash_table->nused--;

	if (hash_table->nused < hash_table->nentries / SHRINK_FACTOR)
	    (void) resize_table(hash_table);
    }

/*
    for (i = 0; i < hash_table->nentries; i++)
    {
        if (hash_table->entries[i].used)
            printf("remove_hash_table2: i = %d, key = %d, data = %d\n",
                        i, (int) (hash_table->entries[i].key),
                        (int) (hash_table->entries[i].data));
    }
*/

    return data;
}

Bool contains_hash_table(Hash_table hash_table, Hash_key key)
{
    Hash_value hash = (hash_table->hash_func)(key);
    Hash_entry *entry = table_entry(hash_table, key, hash);

    return entry && entry->used;
}

Bool find_hash_table(Hash_table hash_table, Hash_key key, Hash_data *data)
{
    Hash_value hash = (hash_table->hash_func)(key);
    Hash_entry *entry = table_entry(hash_table, key, hash);

    if (data && entry && entry->used)
	*data = entry->data;

    return entry && entry->used;
}

void clear_hash_table(Hash_table hash_table,
			Hash_clear_func func, Generic_ptr user_data)
{
    int i;
    Hash_entry *entry;

    for (i = 0; i < hash_table->nentries; i++)
    {
	entry = &(hash_table->entries[i]);

        if (entry->used)
	{
	    entry->used = CCPN_FALSE;
	    if (func)
		(func)(entry->key, entry->data, user_data);
	}
    }

    hash_table->nused = 0;
    (void) resize_table(hash_table);
}

CcpnStatus iterate_hash_table(Hash_table hash_table,
		Hash_iter_func func, Generic_ptr user_data, CcpnString error_msg)
{
    int i;
    Hash_entry *entry;

    for (i = 0; i < hash_table->nentries; i++)
    {
	entry = &(hash_table->entries[i]);

        if (entry->used)
	    CHECK_STATUS((func)(entry->key, entry->data, user_data, error_msg));
    }

    return CCPN_OK;    
}

Hash_value string_hash(int len, char *value)
{
    int i, j, k;
    Hash_value hash = 0, h, v;

    for (i = 0; i < len; i += 4)
    {
	h = 0;
	for (j = i; j < MIN(len, j+4); j++)
	{
	    k = 8 * (j % 4);
	    v = (Hash_value) value[i];
	    h = h | v << k;
	}

	hash = combine_hash(hash, h);
    }

    return hash;
}

Hash_value array_hash(int len, int *values)
{
    int i;
    Hash_value hash = 0, h;
 
    for (i = 0; i < len; i++)
    {
	h = create_hash(values[i]);
	hash = combine_hash(hash, h);
    }

    return hash;
}

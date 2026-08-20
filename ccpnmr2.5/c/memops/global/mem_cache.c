
/*
======================COPYRIGHT/LICENSE START==========================

mem_cache.c: Part of the CcpNmr Analysis program

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
#include "mem_cache.h"

typedef struct Cache_entry
{
    Generic_ptr object;
    int size;
    Cache_delete_func delete_func;
    Generic_ptr delete_data;
    int lock_count;
}   *Cache_entry;

static Cache_entry new_cache_entry(Generic_ptr object, int size,
			Cache_delete_func delete_func, Generic_ptr delete_data)
{
    Cache_entry cache_entry;

    MALLOC_NEW(cache_entry, struct Cache_entry, 1);

    cache_entry->object = object;
    cache_entry->size = size;
    cache_entry->delete_func = delete_func;
    cache_entry->delete_data = delete_data;
    cache_entry->lock_count = 0;

    return cache_entry;
}

static void delete_cache_entry(Cache_entry cache_entry)
{
    if (cache_entry->delete_func)
	(cache_entry->delete_func)(cache_entry->object, cache_entry->delete_data);

    FREE(cache_entry, struct Cache_entry);
}

static void cache_entry_delete(Generic_ptr entry)
{
    delete_cache_entry((Cache_entry) entry);
}

static void lock_cache_entry(Cache_entry cache_entry)
{
    cache_entry->lock_count++;
}

static void unlock_cache_entry(Cache_entry cache_entry)
{
    cache_entry->lock_count--;
}

static Bool is_locked_cache_entry(Cache_entry cache_entry)
{
    if (cache_entry->lock_count > 0)
	return CCPN_TRUE;
    else
	return CCPN_FALSE;
}

static CcpnStatus remove_unlocked(Generic_ptr list_data, Generic_ptr func_data, Bool *done)
{
    Cache_entry entry = (Cache_entry) list_data;
    Mem_cache mem_cache = (Mem_cache) func_data;
    int threshold = (7 * mem_cache->max_size) / 10;

    *done = CCPN_FALSE;

/*
    printf("remove_unlocked: object = %d\n", (int) entry->object);
    printf("remove_unlocked: size = %d\n", entry->size);
    printf("remove_unlocked: lock_count = %d\n", entry->lock_count);
*/

    if (!is_locked_cache_entry(entry))
    {

	CHECK_STATUS(remove_mem_cache(mem_cache, entry->object));

	if (mem_cache->size <= threshold)
	    *done = CCPN_TRUE;
    }

    return CCPN_OK;
}

static void make_room(Mem_cache mem_cache)
{
    Bool forward = CCPN_FALSE;

/*
    printf("make_room: mem_cache->max_size = %d\n", mem_cache->max_size);
    printf("make_room: mem_cache->size = %d\n", mem_cache->size);
*/

    /* should be locked (look at code which calls this function) */
    if (unlock_mutex(mem_cache->mutex) == CCPN_ERROR)
	return;

    traverse_hash_list(mem_cache->list, forward, remove_unlocked, (Generic_ptr) mem_cache);

    /* not much can do if the below fails */
    (void) lock_mutex(mem_cache->mutex);
}

Mem_cache new_mem_cache(int max_size, Equal_func equal_func, Hash_func hash_func)
{
    Mem_cache mem_cache;

    MALLOC_NEW(mem_cache, struct Mem_cache, 1);

    if (!(mem_cache->list = new_hash_list(equal_func, hash_func)))
    {
	FREE(mem_cache, struct Mem_cache);
	return NULL;
    }

    if (!(mem_cache->mutex = new_mutex()))
    {
	delete_mem_cache(mem_cache);
	return NULL;
    }

    mem_cache->max_size = max_size;
    mem_cache->size = 0;

    return mem_cache;
}

void delete_mem_cache(Mem_cache mem_cache)
{
    delete_mutex(mem_cache->mutex);
    delete_hash_list(mem_cache->list, cache_entry_delete);
    FREE(mem_cache, struct Mem_cache);
}

void clear_mem_cache(Mem_cache mem_cache)
{
    clear_hash_list(mem_cache->list, cache_entry_delete);
}

CcpnStatus add_mem_cache(Mem_cache mem_cache, Generic_ptr object, int size,
			Cache_delete_func delete_func, Generic_ptr delete_data)
{
    Cache_entry entry = new_cache_entry(object, size, delete_func, delete_data);

    if (!entry)
	return CCPN_ERROR;

    if (lock_mutex(mem_cache->mutex) == CCPN_ERROR)
    {
	delete_cache_entry(entry);
	return CCPN_ERROR;
    }

    if (insert_hash_list(mem_cache->list, object, (Generic_ptr) entry) == CCPN_ERROR)
    {
	delete_cache_entry(entry);
	(void) unlock_mutex(mem_cache->mutex);
	return CCPN_ERROR;
    }

    mem_cache->size += size;

/*
    printf("add_mem_cache: mem_cache->max_size = %d\n", mem_cache->max_size);
    printf("add_mem_cache: mem_cache->size = %d\n", mem_cache->size);
*/

    lock_cache_entry(entry);

    if (mem_cache->size > mem_cache->max_size)
	make_room(mem_cache);

    CHECK_STATUS(unlock_mutex(mem_cache->mutex));

    return CCPN_OK;
}

CcpnStatus remove_mem_cache(Mem_cache mem_cache, Generic_ptr object)
{
    Cache_entry entry;

    CHECK_STATUS(lock_mutex(mem_cache->mutex));

    entry = (Cache_entry) remove_key_hash_list(mem_cache->list, object);

    if (entry && !is_locked_cache_entry(entry))
    {
	mem_cache->size -= entry->size;
	delete_cache_entry(entry);
    }

    CHECK_STATUS(unlock_mutex(mem_cache->mutex));

    return CCPN_OK;
}

CcpnStatus lock_mem_cache(Mem_cache mem_cache, Generic_ptr object)
{
    Cache_entry entry;

    CHECK_STATUS(lock_mutex(mem_cache->mutex));

    entry = (Cache_entry) promote_hash_list(mem_cache->list, object);

    if (entry)
	lock_cache_entry(entry);

    CHECK_STATUS(unlock_mutex(mem_cache->mutex));

    return CCPN_OK;
}

CcpnStatus unlock_mem_cache(Mem_cache mem_cache, Generic_ptr object)
{
    Cache_entry entry;

    CHECK_STATUS(lock_mutex(mem_cache->mutex));

    entry = (Cache_entry) find_key_hash_list(mem_cache->list, object);

    if (entry)
	unlock_cache_entry(entry);

    CHECK_STATUS(unlock_mutex(mem_cache->mutex));

    return CCPN_OK;
}

CcpnStatus resize_mem_cache(Mem_cache mem_cache, int max_size)
{
    CHECK_STATUS(lock_mutex(mem_cache->mutex));

    mem_cache->max_size = max_size;

    if (mem_cache->size > max_size);
	make_room(mem_cache);

    CHECK_STATUS(unlock_mutex(mem_cache->mutex));

    return CCPN_OK;
}

static CcpnStatus check_func(Generic_ptr list_data, Generic_ptr func_data, Bool *done)
{
    Cache_entry entry = (Cache_entry) list_data;

    printf("mem_cache: check_func: %x: %d: %d\n", (unsigned int) entry->object, entry->size, entry->lock_count);

    return CCPN_OK;
}

CcpnStatus check_mem_cache(Mem_cache mem_cache, CcpnString error_msg)
{
    printf("check_mem_cache\n");

    CHECK_STATUS(traverse_hash_list(mem_cache->list, CCPN_TRUE, check_func, NULL));

    return CCPN_OK;
}


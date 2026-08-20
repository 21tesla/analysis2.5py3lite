
/*
======================COPYRIGHT/LICENSE START==========================

hash_list.c: Part of the CcpNmr Analysis program

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
/* The following code is a slightly modified version of that found in
   "Data structures, an advanced approach using C" by J. Esakov and
   T. Weiss (Prentice Hall, 1989)                                     */

#include "hash_list.h"

typedef struct Hash_node
{
    Generic_ptr data;
    Hash_key key;
    struct Hash_node *next;
    struct Hash_node *prev;
}   *Hash_node;

static Hash_node new_hash_node(Generic_ptr data, Hash_key key)
{
    Hash_node hash_node = (Hash_node) NULL;

    MALLOC_NEW(hash_node, struct Hash_node, 1);

    hash_node->data = data;
    hash_node->key = key;
    hash_node->next = hash_node->prev = (Hash_node) NULL;

    return hash_node;
}
 
static void delete_hash_node(Hash_node hash_node)
{
    FREE(hash_node, struct Hash_node);
}
 
static void delete_hash_nodes(Hash_node hash_node, Data_func func)
{
    Hash_node node;

    for (node = hash_node; node; node = node->next)
    {
        if (func != NULL)  (func)(node->data);
        delete_hash_node(node);
    }
}
 
Hash_list new_hash_list(Equal_func equal_func, Hash_func hash_func)
{
    Hash_list hash_list = (Hash_list) NULL;

    MALLOC_NEW(hash_list, struct Hash_list, 1);

    hash_list->table = new_hash_table(equal_func, hash_func);

    if (!hash_list->table)
    {
	FREE(hash_list, struct Hash_list);
	return NULL;
    }

    hash_list->first = hash_list->last = (Hash_node) NULL;

    return hash_list;
}
 
void delete_hash_list(Hash_list hash_list, Data_func func)
{
    if (hash_list->first)
        delete_hash_nodes(hash_list->first, func);
    delete_hash_table(hash_list->table);
    FREE(hash_list, struct Hash_list);
}

void clear_hash_list(Hash_list hash_list, Data_func func)
{
    clear_hash_table(hash_list->table, (Hash_clear_func) NULL, (Generic_ptr) NULL);
    if (hash_list->first)
    {
        delete_hash_nodes(hash_list->first, func);
        hash_list->first = hash_list->last = (Hash_node) NULL;
    }
}

Bool empty_hash_list(Hash_list hash_list)
{
    return hash_list->first ? CCPN_FALSE : CCPN_TRUE;
}
 
static void set_first_hash_list(Hash_list hash_list, Hash_node hash_node)
{
    if (hash_list->first)
	hash_list->first->prev = hash_node;

    hash_node->next = hash_list->first;
    hash_list->first = hash_node;

    if (!hash_list->last)
	hash_list->last = hash_node;

    hash_node->prev = (Hash_node) NULL;
}

CcpnStatus insert_hash_list(Hash_list hash_list, Hash_key key, Hash_data data)
{
    Hash_node hash_node;
 
    hash_node = new_hash_node(data, key);
 
    if (!hash_node)
	return CCPN_ERROR;

    if (insert_hash_table(hash_list->table, key, (Hash_data) hash_node) == CCPN_ERROR)
    {
	delete_hash_node(hash_node);
	return CCPN_ERROR;
    }

    set_first_hash_list(hash_list, hash_node);

    return  CCPN_OK;
}

static void set_last_hash_list(Hash_list hash_list, Hash_node hash_node)
{
    if (hash_list->last)
	hash_list->last->next = hash_node;

    hash_node->prev = hash_list->last;
    hash_list->last = hash_node;

    if (!hash_list->first)
	hash_list->first = hash_node;

    hash_node->next = (Hash_node) NULL;
}

CcpnStatus append_hash_list(Hash_list hash_list, Hash_key key, Hash_data data)
{
    Hash_node hash_node;
 
    hash_node = new_hash_node(data, key);
 
    if (!hash_node)
	return CCPN_ERROR;

    if (insert_hash_table(hash_list->table, key, (Hash_data) hash_node) == CCPN_ERROR)
    {
	delete_hash_node(hash_node);
	return CCPN_ERROR;
    }

    set_last_hash_list(hash_list, hash_node);

    return  CCPN_OK;
}

static Hash_node find_key_hash_node(Hash_list hash_list, Hash_key key)
{
    Hash_data data;
    Bool found = find_hash_table(hash_list->table, key, &data);

    if (found)
	return (Hash_node) data;
    else
	return (Hash_node) NULL;
}

Generic_ptr find_key_hash_list(Hash_list hash_list, Hash_key key)
{
    Hash_node hash_node = find_key_hash_node(hash_list, key);

    if (hash_node)
	return hash_node->data;
    else
	return (Generic_ptr) NULL;
}

static void cut_hash_list(Hash_list hash_list, Hash_node hash_node)
{
    if (hash_node->prev)
	hash_node->prev->next = hash_node->next;
    else
	hash_list->first = hash_node->next;

    if (hash_node->next)
	hash_node->next->prev = hash_node->prev;
    else
	hash_list->last = hash_node->prev;

    hash_node->prev = hash_node->next = (Hash_node) NULL;
}

static Generic_ptr remove_hash_list(Hash_list hash_list, Hash_node hash_node)
{
    Generic_ptr data = hash_node->data;

    cut_hash_list(hash_list, hash_node);
    (void) remove_hash_table(hash_list->table, hash_node->key);
    delete_hash_node(hash_node);

    return data;
}

Generic_ptr remove_key_hash_list(Hash_list hash_list, Hash_key key)
{
    Hash_data data;
    Bool found = find_hash_table(hash_list->table, key, &data);

    if (found)
	return remove_hash_list(hash_list, (Hash_node) data);
    else
	return (Generic_ptr) NULL;
}

Generic_ptr remove_first_hash_list(Hash_list hash_list)
{
    return remove_hash_list(hash_list, hash_list->first);
}

Generic_ptr remove_last_hash_list(Hash_list hash_list)
{
    return remove_hash_list(hash_list, hash_list->last);
}

CcpnStatus traverse_hash_list(Hash_list hash_list, Bool forward,
			Hash_list_func func, Generic_ptr data)
{
    Bool done = CCPN_FALSE;
    Hash_node hash_node, next_hash_node;

    if (forward)
	hash_node = hash_list->first;
    else
	hash_node = hash_list->last;

    while (hash_node)
    {
	if (forward)
	    next_hash_node = hash_node->next;
	else
	    next_hash_node = hash_node->prev;

	CHECK_STATUS((func)(hash_node->data, data, &done));

	if (done)
	    break;

        hash_node = next_hash_node;
    }

    return CCPN_OK;
}
 
Generic_ptr promote_hash_list(Hash_list hash_list, Hash_key key)
{
    Hash_node hash_node = find_key_hash_node(hash_list, key);

    if (hash_node)
    {
	if (hash_list->first != hash_node)
	{
	    cut_hash_list(hash_list, hash_node);
	    set_first_hash_list(hash_list, hash_node);
	}

	return hash_node->data;
    }

    return (Generic_ptr) NULL;
}

Generic_ptr demote_hash_list(Hash_list hash_list, Hash_key key)
{
    Hash_node hash_node = find_key_hash_node(hash_list, key);

    if (hash_node)
    {
	if (hash_list->last != hash_node)
	{
	    cut_hash_list(hash_list, hash_node);
	    set_last_hash_list(hash_list, hash_node);
	}

	return hash_node->data;
    }

    return (Generic_ptr) NULL;
}


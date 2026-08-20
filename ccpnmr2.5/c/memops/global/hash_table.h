
/*
======================COPYRIGHT/LICENSE START==========================

hash_table.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_hash_table
#define _incl_hash_table

#include "macros.h"
#include "types.h"

typedef Generic_ptr Hash_key;
typedef Generic_ptr Hash_data;

struct Hash_entry;

typedef struct Hash_table
{
    int nentries;
    int nused;
    struct Hash_entry *entries;
    Equal_func equal_func;
    Hash_func hash_func;
}   *Hash_table;

typedef void (*Hash_clear_func)(Hash_key key, Hash_data data,
						Generic_ptr user_data);

typedef CcpnStatus (*Hash_iter_func)(Hash_key key, Hash_data data,
				Generic_ptr user_data, CcpnString error_msg);

extern Bool equal_pointers(Hash_key key1, Hash_key key2);
    
extern Hash_value hash_pointers(Hash_key key);

extern Hash_table new_hash_table(Equal_func equal_func, Hash_func hash_func);

extern void delete_hash_table(Hash_table hash_table);

extern CcpnStatus insert_hash_table(Hash_table hash_table, Hash_key key,
							Hash_data data);

extern Hash_data remove_hash_table(Hash_table hash_table, Hash_key key);

extern Bool contains_hash_table(Hash_table hash_table, Hash_key key);

extern Bool find_hash_table(Hash_table hash_table, Hash_key key,
							Hash_data *data);

extern void clear_hash_table(Hash_table hash_table,
			Hash_clear_func func, Generic_ptr user_data);

extern CcpnStatus iterate_hash_table(Hash_table hash_table,
		Hash_iter_func func, Generic_ptr user_data, CcpnString error_msg);

extern Hash_value string_hash(int len, char *value);

extern Hash_value array_hash(int len, int *values);

#endif /* _incl_hash_table */

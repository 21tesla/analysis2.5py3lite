
/*
======================COPYRIGHT/LICENSE START==========================

hash_list.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_hash_list
#define _incl_hash_list

#include "macros.h"
#include "types.h"

#include "hash_table.h"

struct Hash_node;

typedef struct Hash_list
{
    struct Hash_node *first;
    struct Hash_node *last;
    Hash_table table;
}   *Hash_list;

typedef CcpnStatus (Hash_list_func)(Generic_ptr list_data, Generic_ptr func_data, Bool *done);

extern Hash_list new_hash_list(Equal_func equal_func, Hash_func hash_func);

extern void delete_hash_list(Hash_list hash_list, Data_func func);

extern void clear_hash_list(Hash_list hash_list, Data_func func);

extern Bool empty_hash_list(Hash_list hash_list);

extern CcpnStatus insert_hash_list(Hash_list hash_list, Hash_key key, Generic_ptr data);

extern CcpnStatus append_hash_list(Hash_list hash_list, Hash_key key, Generic_ptr data);

extern Generic_ptr find_key_hash_list(Hash_list hash_list, Hash_key key);

extern Generic_ptr remove_key_hash_list(Hash_list hash_list, Hash_key key);

extern Generic_ptr remove_first_hash_list(Hash_list hash_list);

extern Generic_ptr remove_last_hash_list(Hash_list hash_list);

extern CcpnStatus traverse_hash_list(Hash_list hash_list, Bool forward,
			Hash_list_func func, Generic_ptr data);

extern Generic_ptr promote_hash_list(Hash_list hash_list, Hash_key key);

extern Generic_ptr demote_hash_list(Hash_list hash_list, Hash_key key);

#endif /* _incl_hash_list */


/*
======================COPYRIGHT/LICENSE START==========================

mem_cache.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_mem_cache
#define _incl_mem_cache

#include "hash_list.h"
#include "mutex.h"

typedef struct Mem_cache
{
    int max_size;
    int size;
    Hash_list list;
    Mutex mutex;
}   *Mem_cache;

typedef void (*Cache_delete_func)(Generic_ptr object, Generic_ptr data);

/* max_size is in bytes */
extern Mem_cache new_mem_cache(int max_size,
			Equal_func equal_func, Hash_func hash_func);

extern void delete_mem_cache(Mem_cache mem_cache);

extern void clear_mem_cache(Mem_cache mem_cache);

/* automatically locks object */
extern CcpnStatus add_mem_cache(Mem_cache mem_cache, Generic_ptr object, int size,
		Cache_delete_func delete_func, Generic_ptr delete_data);

extern CcpnStatus remove_mem_cache(Mem_cache mem_cache, Generic_ptr object);

/* promotes object to front of list (so last removed when freeing memory) */
extern CcpnStatus lock_mem_cache(Mem_cache mem_cache, Generic_ptr object);

extern CcpnStatus unlock_mem_cache(Mem_cache mem_cache, Generic_ptr object);

extern CcpnStatus resize_mem_cache(Mem_cache mem_cache, int max_size);

extern CcpnStatus check_mem_cache(Mem_cache mem_cache, CcpnString error_msg);

#endif /* _incl_mem_cache */


/*
======================COPYRIGHT/LICENSE START==========================

list.h: Part of the CcpNmr Analysis program

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
#ifndef _incl_list
#define _incl_list

#include "macros.h"
#include "types.h"

#define  DATA(L)  ((L)->data_ptr)
#define  NEXT(L)  ((L)->next)

#define  CHECK_INSERT_LIST(list, data_ptr, type) \
	 {   if (insert_list(&list, (Generic_ptr) data_ptr) == CCPN_ERROR) \
	     {   FREE(data_ptr, type); \
		 return  CCPN_ERROR;   }   }

typedef struct node
{
    Generic_ptr data_ptr;
    struct node *next;
}   Node, *List;

extern CcpnStatus allocate_node
	(List *p_list, Generic_ptr data_ptr);

extern void free_node
	(List *p_list);

extern CcpnStatus init_list
	(List *p_list);

extern Bool empty_list
	(List list);

extern CcpnStatus insert_list
	(List *p_list, Generic_ptr data_ptr);

extern CcpnStatus append_list
	(List *p_list, Generic_ptr data_ptr);

extern CcpnStatus delete_key
	(List *p_list, Generic_ptr key_ptr);

extern CcpnStatus delete_node
	(List *p_list, List list);

extern CcpnStatus delete_list
	(List *p_list, Generic_ptr *p_data_ptr);

extern CcpnStatus traverse_list
	(List list, Check_func func);

extern List list_iterator
	(List list, List last_return);

extern CcpnStatus find_key
	(List list, Generic_ptr key_ptr, Equal_func cmp, List *p_key_list);

extern void destroy_list
	(List *p_list, Data_func func);

#endif /* _incl_list */


/*
======================COPYRIGHT/LICENSE START==========================

list.c: Part of the CcpNmr Analysis program

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

#include "list.h"

CcpnStatus allocate_node(List *p_list, Generic_ptr data_ptr)
{
    MALLOC(*p_list, Node, 1);
    DATA(*p_list) = data_ptr;
    NEXT(*p_list) = NULL;

    return  CCPN_OK;
}
 
void free_node(List *p_list)
{
    FREE(*p_list, Node);
}
 
CcpnStatus init_list(List *p_list)
{
    *p_list = NULL;
    return  CCPN_OK;
}
 
Bool empty_list(List list)
{
    return  (list == NULL) ? CCPN_TRUE : CCPN_FALSE;
}
 
CcpnStatus insert_list(List *p_list, Generic_ptr data_ptr)
{
    List list;
 
    CHECK_STATUS(allocate_node(&list, data_ptr));
 
    NEXT(list) = *p_list;
    *p_list = list;
    return  CCPN_OK;
}
 
CcpnStatus append_list(List *p_list, Generic_ptr data_ptr)
{
    List list, L;
 
    CHECK_STATUS(allocate_node(&list, data_ptr));
 
    if (empty_list(*p_list) == CCPN_TRUE)
        *p_list = list;
    else
    {
        for (L = *p_list; NEXT(L) != NULL; L = NEXT(L))  ;

        NEXT(L) = list;
    }

    return  CCPN_OK;
}
 
CcpnStatus delete_key(List *p_list, Generic_ptr key_ptr)
{
    List list;

    if (empty_list(*p_list) == CCPN_TRUE)  return  CCPN_ERROR;
 
    if (DATA(*p_list) == key_ptr)
    {
	list = *p_list;
        *p_list = NEXT(*p_list);
    }
    else
    {
        List L;
 
        for (L = *p_list; (L != NULL) &&
			(DATA(NEXT(L)) != key_ptr); L = NEXT(L))  ;
 
        if (L == NULL)  return  CCPN_ERROR;
	list = NEXT(L);
        NEXT(L) = NEXT(list);
    }

    free_node(&list);
 
    return  CCPN_OK;
}
 
CcpnStatus delete_node(List *p_list, List list)
{
    if (empty_list(*p_list) == CCPN_TRUE)  return  CCPN_ERROR;
 
    if (*p_list == list)
    {
        *p_list = NEXT(*p_list);
    }
    else
    {
        List L;
 
        for (L = *p_list; (L != NULL) && (NEXT(L) != list); L = NEXT(L))  ;
 
        if (L == NULL)  return  CCPN_ERROR;
        NEXT(L) = NEXT(list);
    }

    free_node(&list);
 
    return  CCPN_OK;
}
 
CcpnStatus delete_list(List *p_list, Generic_ptr *p_data_ptr)
{
    if (empty_list(*p_list) == CCPN_TRUE)  return  CCPN_ERROR;
 
    *p_data_ptr = DATA(*p_list);
    return  delete_node(p_list, *p_list);
}
 
CcpnStatus traverse_list(List list, Check_func func)
{
    if (empty_list(list) == CCPN_TRUE)  return  CCPN_OK;
 
    CHECK_STATUS((*func)(DATA(list)));
 
    return  traverse_list(NEXT(list), func);
}
 
List list_iterator(List list, List last_return)
{
    return  (last_return == NULL)  ?  list  :  NEXT(last_return);
}
 
CcpnStatus find_key(List list, Generic_ptr key_ptr, Equal_func cmp,
							List *p_key_list)
{
    List curr = NULL;
 
    while ((curr = list_iterator(list, curr)) != NULL)
    {
        if ((*cmp)(key_ptr, DATA(curr)) == CCPN_TRUE)
        {
            *p_key_list = curr;
            return  CCPN_OK;
        }
    }
 
    return  CCPN_ERROR;
}
 
void destroy_list(List *p_list, Data_func func)
{
    if (empty_list(*p_list) == CCPN_FALSE)
    {
        destroy_list(&NEXT(*p_list), func);
        if (func != NULL)  (*func)(DATA(*p_list));
/*        if (func != NULL)  (*func)(&DATA(*p_list)); */  /* modified (&)*/
        free_node(p_list);
    }
}

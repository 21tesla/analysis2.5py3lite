/*
======================COPYRIGHT/LICENSE START==========================

sorts.c: Part of the CcpNmr Analysis program

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
#include "sorts.h"

void heap_sort(Generic_ptr *pointer, int n, Bool ascending,
					Compare_func cmp_fnc)
{
    int c, i, j, m, p;
    Generic_ptr ptr;

    if (ascending)
	c = -1;
    else
	c = 1;

    p = n / 2;
    m = n - 1;

    while (m > 0)
    {
	if (p > 0)
	{
	    p--;
	    ptr = pointer[p];
	}
	else
	{
	    ptr = pointer[m];
	    pointer[m] = pointer[0];
	    m--;
        }

	i = p;
	j = 2*p + 1;

	while (j < m)
	{
	    if (cmp_fnc(pointer[j], pointer[j+1]) == c)
		j++;

	    if (cmp_fnc(ptr, pointer[j]) == c)
	    {
		pointer[i] = pointer[j];
		i = j;
                j = 2*i + 1;
	    }
	    else
	    {
            	j = m + 1;
            }
	}

	if ((j == m) && (cmp_fnc(ptr, pointer[m]) == c))
	{
	    pointer[i] = pointer[m];
            i = m;
        }

        pointer[i] = ptr;
    }
}

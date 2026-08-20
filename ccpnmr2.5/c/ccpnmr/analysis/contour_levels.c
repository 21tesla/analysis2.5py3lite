/*
======================COPYRIGHT/LICENSE START==========================

contour_levels.c: Part of the CcpNmr Analysis program

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
#include "contour_levels.h"

Contour_levels new_contour_levels(int nlevels, float *levels)
{
    int i, j;
    Contour_levels contour_levels;

    MALLOC_NEW(contour_levels, struct Contour_levels, 1);

    MALLOC_NEW(contour_levels->levels, float, nlevels);

    contour_levels->nlevels = nlevels;

    /* first negative contours, then positive ones */
    j = 0;
    for (i = 0; i <  nlevels; i++)
    {
        if (levels[i] < 0)
            contour_levels->levels[j++] = levels[i];
    }

    for (i = 0; i <  nlevels; i++)
    {
        if (levels[i] >= 0)
            contour_levels->levels[j++] = levels[i];
    }

    return contour_levels;
}

void delete_contour_levels(Contour_levels contour_levels)
{
    if (contour_levels)
        FREE(contour_levels->levels, float);

    FREE(contour_levels, struct Contour_levels);
}

Contour_levels copy_contour_levels(Contour_levels contour_levels)
{
    return new_contour_levels(contour_levels->nlevels, contour_levels->levels);
}

Bool equal_contour_levels(Contour_levels contour_levels1,
                          Contour_levels contour_levels2)
{
    int i;

    if (contour_levels1->nlevels != contour_levels2->nlevels)
        return CCPN_FALSE;

    for (i = 0; i < contour_levels1->nlevels; i++)
    {
        if (contour_levels1->levels[i] != contour_levels2->levels[i])
            return CCPN_FALSE;
    }

    return CCPN_TRUE;
}

Bool have_neg_contour_levels(Contour_levels contour_levels)
{
    int i;

    for (i = 0; i < contour_levels->nlevels; i++)
    {
        if (contour_levels->levels[i] < 0)
            return CCPN_TRUE;
    }

    return CCPN_FALSE;
}

Bool have_pos_contour_levels(Contour_levels contour_levels)
{
    int i;

    for (i = 0; i < contour_levels->nlevels; i++)
    {
        if (contour_levels->levels[i] >= 0)
            return CCPN_TRUE;
    }

    return CCPN_FALSE;
}


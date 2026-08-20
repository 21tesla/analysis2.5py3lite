
/*
======================COPYRIGHT/LICENSE START==========================

clipping.c: Part of the CcpNmr Analysis program

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

#include "clipping.h"

typedef enum { POLYLINE_ALL_DRAWN, POLYLINE_NONE_DRAWN, POLYLINE_SOME_DRAWN } PolylineDrawn;

static PolylineDrawn polyline_drawn(Poly_line polyline,
				float xmin, float ymin, float xmax, float ymax)
{
    int i, n = polyline->nvertices;
    Point2f *v = polyline->vertices;
    Bool xlow = CCPN_TRUE, xhigh = CCPN_TRUE, xmid = CCPN_TRUE;
    Bool ylow = CCPN_TRUE, yhigh = CCPN_TRUE, ymid = CCPN_TRUE;

    if (n == 0)
	return POLYLINE_NONE_DRAWN;

    for (i = 0; i < n; i++)
    {
	if (v[i].x > xmin)
	    xlow = CCPN_FALSE;

	if (v[i].x < xmax)
	    xhigh = CCPN_FALSE;

	if ((v[i].x < xmin) || (v[i].x > xmax))
	    xmid = CCPN_FALSE;

	if (v[i].y > ymin)
	    ylow = CCPN_FALSE;

	if (v[i].y < ymax)
	    yhigh = CCPN_FALSE;

	if ((v[i].y < ymin) || (v[i].y > ymax))
	    ymid = CCPN_FALSE;
    }

    if (xmid && ymid)
	return POLYLINE_ALL_DRAWN;

    if (xhigh || xlow || yhigh || ylow)
	return POLYLINE_NONE_DRAWN;

    /* it's possible none might be drawn but have to investigate further */

    return POLYLINE_SOME_DRAWN;
}

void draw_clipped_line(float x0, float y0, float x1, float y1,
		Drawing_funcs *drawing_funcs, Generic_ptr data,
		float xmin, float ymin, float xmax, float ymax)
{
    Bool xclipped = CCPN_TRUE, yclipped = CCPN_TRUE;
    int end_ind = 0, ind0, ind1, min_ind, max_ind, i;
    float xends[4], yends[4], c0, d0, c1, d1, min_x, max_x;

    if (x0 < x1)
    {
	if ((x1 < xmin) || (x0 > xmax))
	    return;

	if ((xmin < x0) & (x1 < xmax))
	    xclipped = CCPN_FALSE;
    }
    else
    {
	if ((x0 < xmin) || (x1 > xmax))
	    return;

	if ((xmin < x1) & (x0 < xmax))
	    xclipped = CCPN_FALSE;
    }

    if (y0 < y1)
    {
	if ((y1 < ymin) || (y0 > ymax))
	    return;

	if ((ymin < y0) & (y1 < ymax))
	    yclipped = CCPN_FALSE;
    }
    else
    {
	if ((y0 < ymin) || (y1 > ymax))
	    return;

	if ((ymin < y1) & (y0 < ymax))
	    yclipped = CCPN_FALSE;
    }

    if (xclipped || yclipped)
    {
	d0 = (xmin-x1)*y0/(x0-x1) + (x0-xmin)*y1/(x0-x1);
	d1 = (xmax-x1)*y0/(x0-x1) + (x0-xmax)*y1/(x0-x1);

	if ((ymin < d0) && (d0 < ymax))
	{
	    xends[end_ind] = xmin;
	    yends[end_ind] = d0;
	    end_ind++;
	}

	if ((ymin < d1) && (d1 < ymax))
	{
	    xends[end_ind] = xmax;
	    yends[end_ind] = d1;
	    end_ind++;
	}

	c0 = (ymin-y1)*x0/(y0-y1) + (y0-ymin)*x1/(y0-y1);
	c1 = (ymax-y1)*x0/(y0-y1) + (y0-ymax)*x1/(y0-y1);

	if ((xmin < c0) && (c0 < xmax))
	{
	    xends[end_ind] = c0;
	    yends[end_ind] = ymin;
	    end_ind++;
	}

	if ((xmin < c1) && (c1 < xmax))
	{
	    xends[end_ind] = c1;
	    yends[end_ind] = ymax;
	    end_ind++;
	}

	if (end_ind != 2)
	    return;

	xends[2] = x0;
	yends[2] = y0;
	xends[3] = x1;
	yends[3] = y1;

	/* below finds middle two indices when sorting on (x, y) */
	/* the middle two indices define the coordinates one needs */
        min_ind = max_ind = 0;
        min_x = max_x = xends[0];
        for (i = 1; i < 4; i++)
	{
	    if (xends[i] < min_x)
	    {
		min_x = xends[i];
		min_ind = i;
	    }
	    else if ((xends[i] == min_x) && (yends[i] < yends[min_ind]))
	    {
		min_x = xends[i];
		min_ind = i;
	    }

	    if (xends[i] > max_x)
	    {
		max_x = xends[i];
		max_ind = i;
	    }
	    else if ((xends[i] == max_x) && (yends[i] > yends[min_ind]))
	    {
		max_x = xends[i];
		max_ind = i;
	    }
	}

        ind0 = ind1 = -1;
        for (i = 0; i < 4; i++)
	{
	    if ((i != min_ind) && (i != max_ind))
	    {
		if (ind0 < 0)
		{
		    ind0 = i;
		}
		else
		{
		    ind1 = i;
		    break;
		}
	    }
	}

	x0 = xends[ind0];
	y0 = yends[ind0];
	x1 = xends[ind1];
	y1 = yends[ind1];
    }

    (drawing_funcs->draw_line)(data, x0, y0, x1, y1);
}

void draw_clipped_polyline(Poly_line polyline, Drawing_funcs *drawing_funcs,
		Generic_ptr data, float xmin, float ymin, float xmax, float ymax)
{
    int i, j, n;
    Point2f *v;
    PolylineDrawn drawn = polyline_drawn(polyline, xmin, ymin, xmax, ymax);

    if (drawn == POLYLINE_ALL_DRAWN)
    {
	(drawing_funcs->draw_polyline)(data, polyline);
    }
    else if (drawn == POLYLINE_SOME_DRAWN)
    {
	n = polyline->nvertices - 1;
	v = polyline->vertices;
	for (i = 0; i < n; i++)
	{
	    j = i + 1;
	    (drawing_funcs->draw_clipped_line)(data, v[i].x, v[i].y, v[j].x, v[j].y);
	}

	if (polyline->closed)
	    (drawing_funcs->draw_clipped_line)(data, v[n].x, v[n].y, v[0].x, v[0].y);
    }
}


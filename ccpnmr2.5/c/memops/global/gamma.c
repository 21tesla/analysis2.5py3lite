/*
======================COPYRIGHT/LICENSE START==========================

gamma.c: Part of the CcpNmr Analysis program

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
#include "gamma.h"

#define ITMAX  100
#define EPS  (3.0e-7)

static double gamma_continued_fraction(double a, double x)
{
    int i;
    double g, gold, a0, a1, b0, b1, f, d;

    g = log_gamma(a);

    gold = 0;
    a0 = 1;
    a1 = x;
    b0 = 0;
    b1 = 1;
    f = 1;

    for (i = 0; i < ITMAX; i++)
    {
	d = i - a;
	a0 = (a1 + d*a0) * f;
	b0 = (b1 + d*a0) * f;

	d = i * f;
	a1 = x*a0 + d*a1;
	b1 = x*b0 + d*b1;

	if (a1 != 0)
	{
	    f = 1.0 / a1;
	    d = b1 * f;

	    if (ABS((d-gold)/d) < EPS)
		return  (d * exp(-x + a*log(x) - g));
	}
    }

    return 1.0; /* arbitrary, really should be an error */
}

static double gamma_series(double a, double x)
{
    int i;
    double g, b, d, s;

    g = log_gamma(a);

    if (x <= 0)
	return 0.0; /* arbitrary, really should be an error */

    b = a;
    d = s = 1 / a;

    for (i = 0; i < ITMAX; i++)
    {
	b += 1;
	d *= x / b;
	s += d;

	if (ABS(d) < (ABS(s) * EPS))
	    return  (s * exp(-x + a*log(x) - g));
    }

    return 0.0; /* arbitrary, really should be an error */
}

double log_gamma(double x)
{
    static double coeff[6] = { 76.18009173, -86.50532033, 24.01409822,
			-1.231739516, 0.120858003e-2, -0.536382e-5 };
    static double c = 2.50662827465;
    int i;
    double t, s;

    x -= 1;
    t = x + 5.5;
    t = (x + 0.5) * log(t) - t;
    s = 1;

    for (i = 0; i < 6; i++)
    {
	x += 1;
	s += coeff[i] / x;
    }

    return  (t + log(c*s));
}

double complement_incomplete_gamma(double a, double x)
{
    if ((a < 0) || (x < 0))
	return  0; /* arbitrary, really should be an error */

    if (x < (a+1))
	return  (1 - gamma_series(a, x));
    else
	return  gamma_continued_fraction(a, x);
}

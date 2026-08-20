/*
======================COPYRIGHT/LICENSE START==========================

line_fit.c: Part of the CcpNmr Analysis program

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
#include "line_fit.h"

#include "gamma.h"

#define  STT_EPS  1.0e-10

CcpnStatus line_fit(int n, float *x, float *y, float *sigma, float *yfit,
		float *p_a, float *p_b, float *p_std_a, float *p_std_b,
		float *corr_ab, float *goodness, CcpnString error_msg)
{
    int i;
    float *t = yfit, a, b, d, s, sx, sy, stt, sig2, std_a, std_b, chi2;

    if (n < 2)
	RETURN_ERROR_MSG("number of parameters < 2")

    if (sigma)
    {
	s = sx = sy = 0;
	for (i = 0; i < n; i++)
	{
	    sig2 = sigma[i] * sigma[i];
	    s += 1.0 / sig2;
	    sx += x[i] / sig2;
	    sy += y[i] / sig2;
	}
    } 
    else
    {
	s = n;
	sx = sy = 0;

	for (i = 0; i < n; i++)
	{
	    sx += x[i];
	    sy += y[i];
	}
    }
/*
printf("sx = %f, sy = %f\n", sx, sy);
*/

    b = 0;

    for (i = 0; i < n; i++)
    {
	t[i] = x[i] - sx / s;

	if (sigma)
	    t[i] /= sigma[i];

	d = t[i] * y[i];

	if (sigma)
	    d /= sigma[i];

	b += d;
    }

    INNER_PRODUCT(stt, t, t, n);

    if (stt < STT_EPS)
	RETURN_ERROR_MSG("x values all the same (it seems)");

    b /= stt;
    a = (sy - b * sx) / s;

    chi2 = 0;
    for (i = 0; i < n; i++)
    {
	yfit[i] = a + b * x[i];
	d = y[i] - yfit[i];

	if (sigma)
	    d /= sigma[i];

	chi2 += d * d;
    }

    std_a = (1.0 + sx*sx/(s*stt)) / s;
    std_b = 1.0 / stt;

    if (!sigma && (n > 2))
    {
	d = chi2 / (n - 2);
	std_a *= d;
	std_b *= d;
    }

    std_a = (float) sqrt((double) std_a);
    std_b = (float) sqrt((double) std_b);

    *p_a = a;
    *p_b = b;
    *p_std_a = std_a;
    *p_std_b = std_b;

/*
printf("yfit[0] = %f, sx = %f, sigma = %d\n", yfit[0], sx, sigma);
*/
    *corr_ab = - sx / (s * stt * std_a * std_b);

    if (sigma)
	*goodness = complement_incomplete_gamma((double) (n-2)/2,
						(double) chi2/2);
    else
	*goodness = chi2;

    return  CCPN_OK;
}

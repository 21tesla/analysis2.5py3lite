/*
======================COPYRIGHT/LICENSE START==========================

method.c: Part of the CcpNmr Analysis program

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
#include "method.h"

#define  SMALL_VALUE  (1.0e-4)

/* u, v, w are three y values in a row for points x-1, x and x+1 */
/* returns the offset of the center from x */
/* implicitly assumes that if v > 0 then v >= u,w and if v < 0 then v <= u,w */
/* b is used for volume calculation only */
static float fit_center3(float u, float v, float w, float *b, Bool use_log)
{
    float c, d;
    Bool is_positive;

    if (v > 0)
        is_positive = CCPN_TRUE;
    else
        is_positive = CCPN_FALSE;

    if (use_log)
    {
        if (!is_positive)
        {
            u = -u;
            v = -v;
            w = -w;
        }

        if ((u > 0) && (w > 0))
        {
            u = (float) log((double) u);
            v = (float) log((double) v);
            w = (float) log((double) w);
        }
        /* otherwise just do non-log fit */
    }

    d = 0.5 * ABS(2*v - u - w);

    if (is_positive)
        *b = d;
    else
        *b = -d;

    if (d > SMALL_VALUE)
    {
        c = 0.25 * ABS(w - u) / d;

        if (is_positive)
        {
            if (w < u)
                c = -c;
        }
        else
        {
            if (w > u)
                c = -c;
        }

        c = MAX(c, -0.499);
        c = MIN(c, 0.499);
    }
    else
    {
        c = 0;
    }

    return c;
}

#define  FUDGE_MULT  0.01

float fit_volume_gaussian3_method(int ndim, float y,
                                  float *ym, float *yp, Bool *dim_done)
{
    int i;
    /*
        float volume, a, b, c;
    */
    /*
        float volume, b, c;
    */
    float volume = 1, b, c;
    float vm, v, vp;

    /* change to log way of calculating volume to try and exclude infinities */
    /* 13 Aug 2004: do not know what the heck the below was trying to do
       so change back to not using logs, only problem here is overall scale */
    for (i = 0; i < ndim; i++)
    {
        if (!dim_done[i])
            continue;

        vm = ym[i];
        v = y;
        vp = yp[i];
        /*
        printf("  i = %d, ym = %3.2f, y = %3.2f, yp = %3.2f\n", i, vm, y, vp);
        */

        if (v > 0)
        {
            if ((vm > v) || (vp > v))
            {
                y = MAX(vm, vp);
                vm = vp = v;
            }
            else
            {
                vm = MAX(vm, FUDGE_MULT*y);
                vp = MAX(vp, FUDGE_MULT*y);
            }
        }
        else
        {
            if ((vm < v) || (vp < v))
            {
                y = MIN(vm, vp);
                vm = vp = v;
            }
            else
            {
                vm = MIN(vm, FUDGE_MULT*y);
                vp = MIN(vp, FUDGE_MULT*y);
            }
        }

        /*
        printf("  i = %d, vm = %3.2f, y = %3.2f, vp = %3.2f\n", i, vm, y, vp);
        */
        c = fit_center3(vm, y, vp, &b, CCPN_TRUE);
        /*
        	c = fit_center3(ym[i], y, yp[i], &b, CCPN_TRUE);
        */
        /*
        	a = y / (float) exp((double) (- b * c * c));
        	volume *= a;
        */
        /*  13 Aug 2004: don't get what the hell this was doing
                volume += b * c * c;
            think it should be below instead
        */
        /*  below should be fit from y = a exp -b(x-c)^2, with x = -1, 0, 1
            do not know what proper x scale is so overall result is out */
        /*  integral if y is (a sqrt(pi/b)) but just ignore a and use y  */
        volume *= sqrt(PI / (double) ABS(b));
        /*
        printf("  i = %d, b = %3.2f, volume = %3.2f\n", i, b, volume);
        */
    }

    /*
        volume = pow((double) volume, (double) 1.0/ndim);
    */
    /*  below was used before 13 Aug 2004
        volume /= ndim;
        volume += (float) log((double) ABS(y));
        / protect against exp overflow and underflow /
        volume = (float) exp((double) MAX(-21.0, MIN(volume, 21.0)));

        volume *= pow(sqrt(PI), (double) ndim);
    */

    volume *= y;
    /*
    printf("at end y = %3.2f, volume = %3.2f\n", y, volume);
    */

    return volume;
}

void fit_center_parabolic_method(int ndim, float *center, float y,
                                 float *ym, float *yp, Bool *dim_done)
{
    int i;
    float b;

    for (i = 0; i < ndim; i++)
        if (dim_done[i])
            center[i] = fit_center3(ym[i], y, yp[i], &b, CCPN_FALSE);
        else
            center[i] = 0;
}

void fit_center_gaussian3_method(int ndim, float *center, float y,
                                 float *ym, float *yp, Bool *dim_done)
{
    int i;
    float b;

    for (i = 0; i < ndim; i++)
        if (dim_done[i])
            center[i] = fit_center3(ym[i], y, yp[i], &b, CCPN_TRUE);
        else
            center[i] = 0;
}

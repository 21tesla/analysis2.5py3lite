/* distribution.c                            *
 *   to contain all the stuff currently in   *
 *  Gaussian.c                               *
 *  Lorentzian.c                             *
 *  Wolfgang.c                               *
 * Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.
 ********************************************
======================COPYRIGHT/LICENSE START==========================

distribution.c: Part of the Bayes program

Copyright (C) 2003-2010 Daniel O'Donovan (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the author: djo35@bioc.cam.ac.uk
=======================================================================
===========================REFERENCE START=============================
===========================REFERENCE END===============================

*/

#include <stdio.h>
#include <math.h>

#include "distribution.h"

#ifndef PI
#define PI 3.14159265358979323836
#endif 

#ifndef ONE_RT_2PI
#define ONE_RT_2PI 0.3989422804
#endif


/*
 * Gaussian function being used: 
 *
 *  G(x) = h * exp [ -1 / 2 * sigma^2 (x - x0)^2]
 *
 * h = the height of the peak, sigma - width, x0 - centre
*/
double Gaussian(double h, double s, double m, double x)
{
    /* volume: * (1 / s) * ONE_RT_2PI */
    return (h * exp( (-1 / (2 * s * s)) * pow(x - m, 2)));
}

/* Currently this seems to be 4x more costly than the 1d function - ??? */
double Gaussian_2d(double h, double s11, double s22, double m1, double m2, int x1, int x2)
{
    double G1, G2;
    /*G1 = intGaussian(1.0, s11, m1, x1);
    G2 = intGaussian(1.0, s22, m2, x2);*/
    G1 = Gaussian(1.0, s11, m1, x1);
    G2 = Gaussian(1.0, s22, m2, x2);
    return h * G1 * G2;
}

double Gaussian_Nd(int ndim, double h, double *s, double *m, int *x)
{
    int i;
    double G = 1.0;
    for (i = 0; i < ndim; i++)
        G *= Gaussian(1.0, s[2*i], m[2*i], x[i]);
    return h * G;
}


/*
 * Lorentzian function being used
 *
 *    L(x) = h * [ (0.5 * Gamma) / ( (x - x0)^2 + (0.5 * Gamma)^2 )]
 *
 *  h added by myself as (I believe) 0 <= L(x) <= 1
 *  Gamma / Sigma - width characteristic
*/

double Lorentzian(double h, double gamma, double x0, double x)
{
    /* Using gamma2 in numerator normalises height to h at maximum */
    double gamma2 = pow(gamma, 2);
    return  h * ( gamma2 / ( gamma2 + pow((x - x0), 2) ));
}

double Lorentzian_ratio(double h, double gamma, double x0, double x)
{
    double gamma2 = h * gamma;
    return  h * ( gamma2 / ( ( pow((x - x0), 2) + gamma2 )));
}

double Lorentzian_2d(double h, double s11, double s22, double m1, double m2, double x1, double x2)
{
    double L1, L2;
/*
    L1 = Lorentzian_ratio(1, s11, m1, x1); 
    L2 = Lorentzian_ratio(1, s22, m2, x2); 
*/
    L1 = Lorentzian(1., s11, m1, x1); 
    L2 = Lorentzian(1., s22, m2, x2); 
    return h * L1 * L2;
}

double Lorentzian_Nd(int ndim, double h, double *s, double *m, int *x)
{
    int i;
    double L = h;
    for (i = 0; i < ndim; i++)
        L *= Lorentzian(1.0, s[2*i], m[2*i], x[i]);
    return L;
}


/*
 * Wolfgang function being used -=- or maybe 'e sub q' / e_q
 *
 *    L(x) = Gamma^(2 / q) * (Gamma^2 + q(x - x0)^2)^(-1/q)
 *
 *  return h * L(x) 
 *
 *  Gamma^(2 / q) to normalise distribution 0 <= height <= 1
 * 
*/

double Wolfgang(double q, double h, double gamma, double x0, double x)
{
    return h * pow(gamma, 2/q) * pow( pow(gamma, 2) + q * pow(x- x0, 2) , (-1/q));
}

double Wolfgang_2d(double q, double h, double s11, double s22, double m1, double m2, double x1, double x2)
{
    double W1, W2;
    W1 = Wolfgang(q, 1, s11, m1, x1);
    W2 = Wolfgang(q, 1, s22, m2, x2);
    return h * W1 * W2;
}

double Wolfgang_Nd(int ndim, double h, double *s, double *m, int *x, double q)
{
    int i;
    double W = 1.0;
    for (i = 0; i < ndim; i++)
        W *= Wolfgang(q, 1.0, s[2*i], m[2*i], x[i]);
    return h * W;
}

/* ============================= 
======================COPYRIGHT/LICENSE START==========================

distribution.h: Part of the Bayes program

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
/* = Header for distribution.c = */
/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.          */
/* ============================= */
#ifndef _incl_distribution
#define _incl_distribution

extern double Gaussian(double h, double s, double m, double x);
extern double Gaussian_2d(double h, double s11, double s22, double m1, double m2, int x1, int x2);
extern double Gaussian_Nd(int ndim, double h, double *s, double *m, int *x);

extern double Lorentzian(double h, double gamma, double x0, double x);
extern double Lorentzian_2d(double h, double s11, double s22, double m1, double m2, double x1, double x2);
extern double Lorentzian_Nd(int ndim, double h, double *s, double *m, int *x);

extern double Wolfgang(double q, double h, double gamma, double x0, double x);
extern double Wolfgang_2d(double q, double h, double s11, double s22, double m1, double m2, double x1, double x2);
extern double Wolfgang_Nd(int ndim, double h, double *s, double *m, int *x, double q);

#endif /* _incl_distribution */
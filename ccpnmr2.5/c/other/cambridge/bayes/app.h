/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. All rights reserved.
 
======================COPYRIGHT/LICENSE START==========================

app.h: Part of the Bayes program

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
/* header file for bayesapp.c */
#ifndef _incl_app
#define _incl_app

#include "userstr.h"

#ifndef PI
#define PI 3.14159265358979323836
#endif

extern int UserBuild(double*, CommonStr*, ObjectStr*, int);

#endif /* _incl_app */
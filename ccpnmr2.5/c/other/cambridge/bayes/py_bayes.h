/* ============================================ 
======================COPYRIGHT/LICENSE START==========================

py_bayes.h: Part of the Bayes program

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
/* = Shamelessly pinched from py_block_file.h = */
/* Created by Daniel O'Donovan on 2008-09-22.
 * Copyright (c) 2008 University of Cambridge. 
 * All rights reserved.                         */
/* ============================================ */

#ifndef _incl_py_bayes
#define _incl_py_bayes

#include "macros.h"
#include "types.h"

#include "bayes_nmr.h"

typedef struct Py_BayesPeakSeparator
{
    PyObject_HEAD
    Bayes bayes;
}   *Py_Bayes;

extern Bool is_py_bayes(PyObject *obj);

#endif /* _incl_py_bayes */

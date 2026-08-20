/*
======================COPYRIGHT/LICENSE START==========================

midge.h: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Alexei Grishaev, Miguel Llinas, Guillermo Bermejo, Wayne Boucher and Tim Stevens (Carnegie Mellon University and University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
M. Madrid, E. Llinas and M. Llinas (1991).
Model-Independent Refinement of Interproton Distances Generated from
H-1-NMR Overhauser Intensities. 
Journal Of Magnetic Resonance 93: 329-346.

A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/
#ifndef _incl_midge
#define _incl_midge

#include "macros.h"
#include "types.h"

#define  NO_TYPE    0
#define  CH3_TYPE   1
#define  CH2_TYPE   2
#define  HARO_TYPE  3
#define  HN_TYPE    4
#define  KHZ_TYPE   5

typedef struct Midge
{
    int n;
    int *nhs;
    int *types;
    double **copy;
    double **aexp0;
    double **sp;
    double *evalues;
    double *work;
}   *Midge;

extern Midge new_midge(int n, int *nhs, int *types);

extern void delete_midge(Midge midge);

extern CcpnStatus run_midge(Midge midge, double **amat, double **rmat,
	int max_iter, float sf, float tmix, float tcor, float rleak,
	Bool n15_labelled, Bool c13_labelled,
	float *err, CcpnString error_msg);

#endif /* _incl_midge */

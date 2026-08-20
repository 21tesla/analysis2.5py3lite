/*
======================COPYRIGHT/LICENSE START==========================

dynamics.h: Part of the CcpNmr Clouds program

Copyright (C) 2005 Alexander Lemak, Miguel Llinas, Wayne Boucher and Tim Stevens (University of Toronto, Carnegie Mellon University and University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================
*/
#ifndef _incl_dynamics
#define _incl_dynamics

#include "atom_coord_list.h"
#include "dist_constraint_list.h"
#include "dist_force.h"

typedef struct Dynamics
{
    float rp_force_const;
    float beta;
    float rmin;
    float drzap;
    float tref;
    float tau;
    float elapsed_time;
    int nsteps;
    int nprint;
}   *Dynamics;

extern void free_dynamics_memory(void);

extern Dynamics new_dynamics(float rp_force_const, float beta, float rmin,
		float drzap, float tref, float tau, float elapsed_time,
		int nsteps, int nprint);

extern void delete_dynamics(Dynamics dynamics);

extern CcpnStatus run_dynamics(Dynamics dynamics,
		Atom_coord_list atom_coord_list,
		Dist_constraint_list noe_list,
		Dist_force noe_force, CcpnString error_msg);

#endif /* _incl_dynamics */

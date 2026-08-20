/*
======================COPYRIGHT/LICENSE START==========================

random.c: Part of the CcpNmr Analysis program

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
#include "random.h"

#define  RAN_ARRAY_SIZE  123

#define  M1  259200
#define  A1  7141
#define  C1  54773

#define  M2  134456
#define  A2  8121
#define  C2  28411

#define  M3  243000
#define  A3  4561
#define  C3  51349

#define  R1  (1.0 / ((float) M1))
#define  R2  (1.0 / ((float) M2))

static int x1, x2, x3;

static float ran_array[RAN_ARRAY_SIZE];

void set_seed(unsigned int seed)
{
    int i;

    x1 = (C1 + seed) % M1;
    x1 = (A1*x1 + C1) % M1;
    x2 = x1 % M2;
    x1 = (A1*x1 + C1) % M1;
    x3 = x1 % M3;

    for (i = 0; i < RAN_ARRAY_SIZE; i++)
    {
	x1 = (A1*x1 + C1) % M1;
	x2 = (A2*x2 + C2) % M2;

	ran_array[i] = (x1 + x2*R2)*R1;
    }
}

float uniform01()
{
    float r;
    unsigned int ind;
    static Bool initialised = CCPN_FALSE;

    if (!initialised)
    {
	set_seed(123456);
	initialised = CCPN_TRUE;
    }

    x1 = (A1*x1 + C1) % M1;
    x2 = (A2*x2 + C2) % M2;
    x3 = (A3*x3 + C3) % M3;

    ind = (RAN_ARRAY_SIZE * x3) / M3;
    r = ran_array[ind];
    ran_array[ind] = (x1 + x2*R2)*R1;

    return  r;
}

float uniform(float a, float b)
{
    return  (a + (b-a) * uniform01());
}

float normal01()
{
    float x1, x2, r, y;
    static float next_y;
    static Bool flag = CCPN_TRUE;

    if (flag)
    {
	flag = CCPN_FALSE;

	do
	{
	    x1 = uniform01();
	}   while (x1 == 0);

	x2 = TWOPI * uniform01();

	r = sqrt(-2.0 * log(x1));

	y = r * cos(x2);
	next_y = r * sin(x2);
    }
    else
    {
	flag = CCPN_TRUE;

	y = next_y;
    }

    return  y;
}

float normal(float mean, float std_dev)
{
    return  (mean + std_dev * normal01());
}

/*
======================COPYRIGHT/LICENSE START==========================

dynamics.c: Part of the CcpNmr Clouds program

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
#include "dynamics.h"

#include "random.h"

/*  DANGER: 7 Dec 09: assume dist_power = 4
                      and soft_exponent = 1
                           and exponent = 2
    At same time changed calc_avg_distance to calc_avg_distance2
    without changing the commented out pow code to do the right thing
*/

#define  CONSTB  0.83144087
#define  TOKELV  1007.0
#define  TOT0    20.458

#define  EPS     1.0e-8

static int natoms = 0;
static int natoms_alloc = 0;
static int *fixed = NULL;
static double *mass = NULL;
static double *x = NULL;
static double *y = NULL;
static double *z = NULL;
static double *xx = NULL;
static double *yy = NULL;
static double *zz = NULL;
static double *vx = NULL;
static double *vy = NULL;
static double *vz = NULL;
static double *ax = NULL;
static double *ay = NULL;
static double *az = NULL;
static double *fx = NULL;
static double *fy = NULL;
static double *fz = NULL;

static int nrp = 0;
static int nrp_alloc = 0;
static int nrp_incr = 10000;
static int *rp_list0 = NULL;
static int *rp_list1 = NULL;

static void free_atom_memory(void)
{
    FREE(mass, double);
    FREE(x, double);
    FREE(y, double);
    FREE(z, double);
    FREE(fixed, int);
    FREE(xx, double);
    FREE(yy, double);
    FREE(zz, double);
    FREE(vx, double);
    FREE(vy, double);
    FREE(vz, double);
    FREE(ax, double);
    FREE(ay, double);
    FREE(az, double);
    FREE(fx, double);
    FREE(fy, double);
    FREE(fz, double);

    natoms_alloc = 0;
}

static CcpnStatus alloc_atom_memory(CcpnString error_msg)
{
    sprintf(error_msg, "allocating atom memory");

    free_atom_memory();

    MALLOC(mass, double, natoms);
    MALLOC(x, double, natoms);
    MALLOC(y, double, natoms);
    MALLOC(z, double, natoms);
    MALLOC(fixed, int, natoms);
    MALLOC(xx, double, natoms);
    MALLOC(yy, double, natoms);
    MALLOC(zz, double, natoms);
    MALLOC(vx, double, natoms);
    MALLOC(vy, double, natoms);
    MALLOC(vz, double, natoms);
    MALLOC(ax, double, natoms);
    MALLOC(ay, double, natoms);
    MALLOC(az, double, natoms);
    MALLOC(fx, double, natoms);
    MALLOC(fy, double, natoms);
    MALLOC(fz, double, natoms);

    natoms_alloc = natoms;

    return CCPN_OK;
}

static void free_rp_memory(void)
{
    FREE(rp_list0, int);
    FREE(rp_list1, int);

    nrp_alloc = 0;
}

static CcpnStatus alloc_rp_memory(CcpnString error_msg)
{
    int n = nrp_alloc + nrp_incr;

    if (nrp_alloc == 0)
    {
	sprintf(error_msg, "allocating rp memory");
	MALLOC(rp_list0, int, n);
	MALLOC(rp_list1, int, n);
    }
    else
    {
	sprintf(error_msg, "reallocating rp memory");
	REALLOC(rp_list0, int, n);
	REALLOC(rp_list1, int, n);
    }

    nrp_alloc = n;

    return CCPN_OK;
}

void free_dynamics_memory(void)
{
    free_atom_memory();
    free_rp_memory();
}

static CcpnStatus calc_rp_list(double rzap, CcpnString error_msg)
{
    int i, j;
    double dx, dy, dz, d2, rzap2 = rzap * rzap;

    nrp = 0;
    for (i = 0; i < natoms-1; i++)
    {
	for (j = i+1; j < natoms; j++)
	{
	    if ( fixed[i] && fixed[j] ) {
                continue;
            }
            dx = x[i] - x[j];
            if (ABS(dx) > rzap)
		continue;

	    dy = y[i] - y[j];
            if (ABS(dy) > rzap)
		continue;

	    dz = z[i] - z[j];
            if (ABS(dz) > rzap)
		continue;

	    d2 = dx*dx + dy*dy + dz*dz;
	    if (d2 > rzap2)
		continue;

	    if (nrp >= nrp_alloc)
		CHECK_STATUS(alloc_rp_memory(error_msg));

	    rp_list0[nrp] = i;
	    rp_list1[nrp] = j;
	    nrp++;
	}
    }

    return CCPN_OK;
}

static double kinetic_energy(void)
{
    int i;
    double d = 0;

    for (i = 0; i < natoms; i++)
	d += mass[i] * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    return 0.5 * d;
}

static void calc_dist_stats(int nconstraints, Dist_constraint *constraints,
				int *nviol, double *rmsd)
{
    int i, j, k, n = 0;
    double d, dmin, dmax, dx, dy, dz, r, s = 0;
    static double thresh = 0.5;

/**** TBD
    for (i = 0; i < nconstraints; i++)
    {
	j = constraints[i]->atom0;
	k = constraints[i]->atom1;
	if (j == k)
	    continue;

	dmin = constraints[i]->dist_lower;
	dmax = constraints[i]->dist_upper;
	dx = x[j] - x[k];
	dy = y[j] - y[k];
	dz = z[j] - z[k];
	r = sqrt(dx*dx + dy*dy + dz*dz);

	if (r <= dmin-thresh)
	    n++;
        else if (r >= dmax+thresh)
	    n++;

	if (r < dmin)
	    d = dmin - r;
        else if (r > dmax)
	    d = r - dmax;
	else
	    d = 0;

	s += d * d;
    }

    *nviol = n;
    *rmsd = sqrt(s / nconstraints);
****/
}

static double inertial_radius(void)
{
    int i;
    double dx, dy, dz, r = 0, t = 0, xc = 0, yc = 0, zc = 0;

    for (i = 0; i < natoms; i++)
    {
	t += mass[i];
	xc += mass[i] * x[i];
	yc += mass[i] * y[i];
	zc += mass[i] * z[i];
    }

    xc /= t;
    yc /= t;
    zc /= t;

    for (i = 0; i < natoms; i++)
    {
	dx = x[i] - xc;
	dy = y[i] - yc;
	dz = z[i] - zc;
	r += mass[i] * (dx*dx + dy*dy + dz*dz);
    }

    r = sqrt(r / t);

    return r;
}

static void update_all(double tref, double tau, double beta)
{
    int i;
    double ek, r, rtau, temp;

    rtau = 0.5 * tau * tau;
    ek = kinetic_energy();
    temp = TOKELV * ek / (3 * natoms);
    temp = MAX(temp, 0.001);
    r = beta * (tref/temp-1.0);

    for (i = 0; i < natoms; i++)
    {
        if (! fixed[i]) {
	    ax[i] = fx[i] / mass[i] + r * vx[i];
	    ay[i] = fy[i] / mass[i] + r * vy[i];
	    az[i] = fz[i] / mass[i] + r * vz[i];
	    x[i] += tau * vx[i] + rtau * ax[i];
	    y[i] += tau * vy[i] + rtau * ay[i];
	    z[i] += tau * vz[i] + rtau * az[i];
	    vx[i] += tau * ax[i];
	    vy[i] += tau * ay[i];
	    vz[i] += tau * az[i];
        }
    }
}

static void update_v(double tref, double tau, double beta)
{
    int i;
    double ek, r, temp;

    ek = kinetic_energy();
    temp = TOKELV * ek / (3 * natoms);
    temp = MAX(temp, 0.001);
    r = beta * (tref/temp-1.0);

    for (i = 0; i < natoms; i++)
    {
	if (!fixed[i]) {
	    vx[i] += 0.5 * tau * (fx[i] / mass[i] + r * vx[i] - ax[i]);
	    vy[i] += 0.5 * tau * (fy[i] / mass[i] + r * vy[i] - ay[i]);
	    vz[i] += 0.5 * tau * (fz[i] / mass[i] + r * vz[i] - az[i]);
	}
    }
}

static double calc_rp_force(double force_const, double rmin)
{
    int i, j, k;
    double dx, dy, dz, d2, dr, rjk, force = 0;
    double rmin2 = rmin * rmin;

    if (force_const == 0)
      return force;

/*
    printf("*****force constant = %f, rmin = %f\n", force_const, rmin);
*/

    for (i = 0; i < nrp; i++)
    {
	j = rp_list0[i];
	k = rp_list1[i];

	dx = x[k] - x[j];
        if (ABS(dx) > rmin)
	    continue;

	dy = y[k] - y[j];
        if (ABS(dy) > rmin)
	    continue;

	dz = z[k] - z[j];
        if (ABS(dz) > rmin)
	    continue;

	d2 = dx*dx + dy*dy + dz*dz;
	if (d2 > rmin2)
	    continue;

	dr = rmin2 - d2;
	force += force_const * dr * dr;

	rjk = 4 * force_const * dr;
	dx *= rjk;
	dy *= rjk;
	dz *= rjk;

	fx[j] -= dx;
	fx[k] += dx;

	fy[j] -= dy;
	fy[k] += dy;

	fz[j] -= dz;
	fz[k] += dz;
    }

    return force;
}

static double power(double x, double y)
{
    if (y == 0)
	return 1;
    else if (y == 1)
	return x;
    else if (y == 2)
	return x*x;
    else
	return pow(x, y);
}

static double calc_avg_distance2(Dist_constraint constraint, float dist_power)
{
    int i, j, k, natom_pairs = constraint->natom_pairs;
/*
    double s = 0, p = -0.5*dist_power, dx, dy, dz, r2;
*/
    double s = 0, dx, dy, dz, r2;

    if (natom_pairs == 1)
    {
	j = constraint->atoms0[0];
	k = constraint->atoms1[0];
	if (j != k)
	{
	    dx = x[j] - x[k];
	    dy = y[j] - y[k];
	    dz = z[j] - z[k];
	    s = dx*dx + dy*dy + dz*dz;
	}
    }
    else
    {
        for (i = 0; i < constraint->natom_pairs; i++)
        {
	    j = constraint->atoms0[i];
	    k = constraint->atoms1[i];
	    if (j != k)
	    {
	        dx = x[j] - x[k];
	        dy = y[j] - y[k];
	        dz = z[j] - z[k];
	        r2 = dx*dx + dy*dy + dz*dz;
/*
	        s += pow(r2, p);
*/
	        s += 1.0 / (r2 * r2);
	    }
        }

        if (s > 0)
/*
	    s = pow(s, -1.0/dist_power);
*/
	    s = 1.0 / sqrt(s);
    }

    return s;
}

static double calc_dist_force(int nconstraints, Dist_constraint *constraints,
		float force_const, float exponent, float soft_exponent,
		float r_switch, float asymptote, float dist_power)
{
    int i, j, k, n, natom_pairs;
    double a, b, da, d, dmin, dmax, dx, dy, dz, exp1, exp2;
    double r, r2, s2, rs1, rs2, rs3, rs4, rjk, ujk, force = 0, t;
    Dist_constraint constraint;

    rs1 = pow(r_switch, exponent);
    rs2 = pow(r_switch, soft_exponent);
    rs3 = pow(r_switch, soft_exponent + 1);
    rs4 = pow(r_switch, exponent + soft_exponent);
    b = (asymptote*rs3 - exponent*rs4) / soft_exponent;
    a = rs1 - asymptote*r_switch - b/rs2;
    exp1 = exponent - 1;
    exp2 = soft_exponent + 1;
/*
    printf("calc_dist_force: %6.2lf %6.2lf\n", force_const, asymptote);
    printf("calc_dist_force: %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf\n",
	rs1, rs2, rs3, rs4, a, b, exp1, exp2);
*/

    for (i = 0; i < nconstraints; i++)
    {
	constraint = constraints[i];

	r2 = calc_avg_distance2(constraint, dist_power);
	if (r2 <= 0)
	    continue;

	dmin = constraint->dist_lower;
	dmax = constraint->dist_upper;

	da = dmax + r_switch;

	if (r2 < dmin*dmin)
	{
	    r2 = MAX(r2, EPS); /* protect against r2 getting near 0 */
	    r = sqrt(r2);
	    d = dmin - r;
	    ujk = force_const * d * d;
/*
	    rjk = force_const * exponent * pow(d, exp1);
*/
	    rjk = force_const * 2 * d;
	}
        else if (r2 > dmax*dmax)
	{
	    r = sqrt(r2);
	    d = r - dmax;

	    if (r <= da)
	    {
		ujk = force_const * d * d;
/*
		rjk = - force_const * exponent * pow(d, exp1);
*/
	        rjk = - force_const * 2 * d;
	    }
	    else
	    {
/*
		ujk = force_const * (a + asymptote*d + b/pow(d, soft_exponent));
		rjk = - force_const * (asymptote - b*soft_exponent/pow(d, exp2));
*/
		ujk = force_const * (a + asymptote*d + b/d);
		rjk = - force_const * (asymptote - b/(d*d));
	    }
	}
	else
	{
	    ujk = rjk = 0;
	    r = 1.0;
	}

	force += ujk;

	natom_pairs = constraints[i]->natom_pairs;
        if (natom_pairs == 1)
	{
	    j = constraint->atoms0[0];
	    k = constraint->atoms1[0];
	    if (j == k)
		continue;
	    t = rjk / r;
	    dx = x[j] - x[k];
	    dy = y[j] - y[k];
	    dz = z[j] - z[k];
	    dx *= t;
	    dy *= t;
	    dz *= t;
	    fx[j] += dx;
	    fx[k] -= dx;
	    fy[j] += dy;
	    fy[k] -= dy;
	    fz[j] += dz;
	    fz[k] -= dz;
	}
	else
	{
            for (n = 0; n < constraints[i]->natom_pairs; n++)
            {
	        j = constraint->atoms0[n];
	        k = constraint->atoms1[n];
	        if (j == k)
		    continue;
	        dx = x[j] - x[k];
	        dy = y[j] - y[k];
	        dz = z[j] - z[k];
                s2 = dx*dx + dy*dy + dz*dz;
/*
	        t = rjk * pow(r, dist_power+1.0) * pow(s2, -dist_power/2.0-1.0);
*/
		t = rjk * r2 * r2 * r / (s2 * s2 * s2);
	        dx *= t;
	        dy *= t;
	        dz *= t;
	        fx[j] += dx;
	        fx[k] -= dx;
	        fy[j] += dy;
	        fy[k] -= dy;
	        fz[j] += dz;
	        fz[k] -= dz;
            }
	}
    }

    return force;
}

Dynamics new_dynamics(float rp_force_const, float beta, float rmin,
		float drzap, float tref, float tau, float elapsed_time,
		int nsteps, int nprint)
{
    Dynamics dynamics;

    MALLOC_NEW(dynamics, struct Dynamics, 1);

    dynamics->rp_force_const = rp_force_const;
    dynamics->beta = beta;
    dynamics->rmin = rmin;
    dynamics->drzap = drzap;
    dynamics->tref = tref;
    dynamics->tau = tau;
    dynamics->elapsed_time = elapsed_time;
    dynamics->nsteps = nsteps;
    dynamics->nprint = nprint;

    return dynamics;
}

void delete_dynamics(Dynamics dynamics)
{
    FREE(dynamics, struct Dynamics);
}

CcpnStatus run_dynamics(Dynamics dynamics,
		Atom_coord_list atom_coord_list,
		Dist_constraint_list noe_list,
		Dist_force noe_force,
                CcpnString error_msg)
{
    int i, j, n, step;
    double d, d2, dist, dx, dy, dz, ek, rmsd, s, tautot0, temp, unoe, urpl;
    double rp_force_const = dynamics->rp_force_const;
    double beta = dynamics->beta;
    double rmin = dynamics->rmin;
    double drzap = dynamics->drzap;
    double rzap = rmin + drzap;
    double tref = dynamics->tref;
    double tau = dynamics->tau;
    double elapsed_time = dynamics->elapsed_time;
    int nsteps = dynamics->nsteps;
    int nprint = dynamics->nprint;
    Atom_coord *atom_coords = atom_coord_list->atom_coords;
    int nnoes = noe_list->ndist_constraints;
    Dist_constraint *noes = noe_list->dist_constraints;

    natoms = atom_coord_list->natom_coords;

    for (i = 0; i < nnoes; i++)
    {
	for (j = 0; j < noes[i]->natom_pairs; j++)
	{
	    n = noes[i]->atoms0[j];
	    if ((n < 0) || (n >= natoms))
	    {
	        sprintf(error_msg,
		    "NOE constraint %d:%d has atom0 = %d, should be in range [0, %d]",
		    i, j, n, natoms-1);
	    }

	    n = noes[i]->atoms1[j];
	    if ((n < 0) || (n >= natoms))
	    {
	        sprintf(error_msg,
		    "NOE constraint %d:%d has atom1 = %d, should be in range [0, %d]",
		    i, j, n, natoms-1);
	    }
	}
    }
/*
    printf("run_dynamics: natoms=%4d, nsteps=%4d, tref=%3.0f, rmin=%3.1f, drzap=%3.1f, beta=%3.1f, tau=%3.1f\n",
	natoms, nsteps, tref, rmin, drzap, beta, tau);
*/

    if (natoms > natoms_alloc)
	CHECK_STATUS(alloc_atom_memory(error_msg));

    for (i = 0; i < natoms; i++)
    {
	mass[i] = atom_coords[i]->mass;
        fixed[i] = atom_coords[i]->isFixed;
	x[i] = atom_coords[i]->x;
	y[i] = atom_coords[i]->y;
	z[i] = atom_coords[i]->z;
/*
	if (i < 3)
	    printf("runA: i = %d, mass = %3.0f, x = %6.1f, y = %6.1f, z = %6.1f\n",
			i, mass[i], x[i], y[i], z[i]);
*/
    }

    beta /= TOT0;
    tautot0 = tau * TOT0;

    d = 0;
    for (i = 0; i < natoms; i++)
    {
        if ( fixed[i] ) {
	  vx[i] = 0.0;
	  vy[i] = 0.0;
	  vz[i] = 0.0;

        } else {

          s = sqrt(CONSTB * tref / mass[i]);
	  vx[i] = s * normal01() / TOT0;
	  vy[i] = s * normal01() / TOT0;
	  vz[i] = s * normal01() / TOT0;

        }
    }

    ek = 3 * natoms * tref / TOKELV;
    d = kinetic_energy();
    d = sqrt(ek / d);
    SCALE_VECTOR(vx, vx, d, natoms);
    SCALE_VECTOR(vy, vy, d, natoms);
    SCALE_VECTOR(vz, vz, d, natoms);

    printf("temperature = %1.0f\n", tref);

    COPY_VECTOR(xx, x, natoms);
    COPY_VECTOR(yy, y, natoms);
    COPY_VECTOR(zz, z, natoms);

/*
    printf("kinetic energy = %lf\n", kinetic_energy());
    for (i = 0; i < natoms; i++)
	printf("%3d %6.1f %6.1f %6.1f %6.3f %6.3f %6.3f\n",
		i, x[i], y[i], z[i], vx[i], vy[i], vz[i]);
*/

    for (step = 0; step < nsteps; step++)
    {
	dist = 0;

	if (step > 0)
	{
	    for (i = 0; i < natoms; i++)
	    {
		dx = x[i] - xx[i];
		dy = y[i] - yy[i];
		dz = z[i] - zz[i];
		d2 = dx*dx + dy*dy + dz*dz;
		dist = MAX(dist, d2);
	    }

	    dist = sqrt(dist);
	}

	if ((step == 0) || (dist > 0.5*drzap))
	{
	    CHECK_STATUS(calc_rp_list(rzap, error_msg));

	    COPY_VECTOR(xx, x, natoms);
	    COPY_VECTOR(yy, y, natoms);
	    COPY_VECTOR(zz, z, natoms);
	}

	if (step == 0)
	{
	    ZERO_VECTOR(fx, natoms);
	    ZERO_VECTOR(fy, natoms);
	    ZERO_VECTOR(fz, natoms);
	    urpl = calc_rp_force(rp_force_const, rmin);
	    unoe = calc_dist_force(nnoes, noes, noe_force->force_const,
			noe_force->exponent, noe_force->soft_exponent,
			noe_force->r_switch, noe_force->asymptote,
			noe_force->dist_power);
	}

	update_all(tref, tautot0, beta);

	ZERO_VECTOR(fx, natoms);
	ZERO_VECTOR(fy, natoms);
	ZERO_VECTOR(fz, natoms);
	urpl = calc_rp_force(rp_force_const, rmin);
	unoe = calc_dist_force(nnoes, noes, noe_force->force_const,
			noe_force->exponent, noe_force->soft_exponent,
			noe_force->r_switch, noe_force->asymptote,
			noe_force->dist_power);

	update_v(tref, tautot0, beta);

	if (step % nprint == 0)
	{
	    ek = kinetic_energy();
	    temp = TOKELV * ek / (3 * natoms);
	    calc_dist_stats(nnoes, noes, &n, &rmsd);
	    d = inertial_radius();

	    printf("%10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %6.2lf %6d %6d\n",
			elapsed_time, temp, urpl, unoe, rmsd, d, n, nrp);
	}

	elapsed_time += tau;
    }

    for (i = 0; i < natoms; i++)
    {

	    atom_coords[i]->x = x[i];
	    atom_coords[i]->y = y[i];
	    atom_coords[i]->z = z[i];

/*
	if (i < 3)
	    printf("runB: i = %d, mass = %3.0f, x = %6.1f, y = %6.1f, z = %6.1f\n",
			i, mass[i], x[i], y[i], z[i]);
*/
    }

    dynamics->elapsed_time = elapsed_time;

    return CCPN_OK;
}

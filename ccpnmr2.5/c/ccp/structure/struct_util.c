
/*
======================COPYRIGHT/LICENSE START==========================

struct_util.c: Part of the CcpNmr Analysis program

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
#include "struct_util.h"

#include "eigenvalue.h"
#include "geometry.h"
#include "linalg.h"

/*
static void print_matrix(CcpnString msg, float **matrix)
{
    int i, j;

    printf("%s:\n", msg);
    for (i = 0; i < COORD_NDIMS; i++)
    {
        for (j = 0; j < COORD_NDIMS; j++)
            printf(" %f", matrix[i][j]);
        printf("\n");
    }

    printf("\n");
}
*/

/* sort into decreasing order */
static void sort_eigenvalues(float *eigenvalues, float **matrix)
{
    int i, j, k;
    float p;

    for (i = 0; i < COORD_NDIMS-1; i++)
    {
        k = i;
        p = eigenvalues[k];
        for (j = i+1; j < COORD_NDIMS; j++)
        {
            if (eigenvalues[j] > p)
            {
                k = j;
                p = eigenvalues[k];
            }
        }

        if (k != i)
        {
            eigenvalues[k] = eigenvalues[i];
            eigenvalues[i] = p;
            for (j = 0; j < COORD_NDIMS; j++)
            {
                p = matrix[j][i];
                matrix[j][i] = matrix[j][k];
                matrix[j][k] = p;
            }
        }
    }
}

/* translate coords to center of mass, using weight if defined */
/* input: coords(atom, dim), weight(atom) */
/* output: coords translated */
void translate_coordinates(int natoms, float **coords, float *weight)
{
    int i, j;
    float t = 0.0, w = 1.0, center[COORD_NDIMS];

    ZERO_VECTOR(center, COORD_NDIMS);

    for (i = 0; i < natoms; i++)
    {
	if (weight)
	    w = weight[i];

	t += w;
        for (j = 0; j < COORD_NDIMS; j++)
	    center[j] += w * coords[i][j];
    }

    for (j = 0; j < COORD_NDIMS; j++)
	center[j] /= t;

    for (i = 0; i < natoms; i++)
        for (j = 0; j < COORD_NDIMS; j++)
	    coords[i][j] -= center[j];
}

/* align coords1 to coords2, after translating to center of mass */
/* input: coords1(atom, dim), coords2(atom, dim), weight(atom) */
/* output: rotation(dim, dim), err */
CcpnStatus align_translate_coordinates(int natoms, float **coords1, float **coords2,
                         float *weight, float **rotation, float *err, CcpnString error_msg)
{
    /* translate to center of mass */

    translate_coordinates(natoms, coords1, weight);
    translate_coordinates(natoms, coords2, weight);

    return align_coordinates(natoms, coords1, coords2, weight, rotation, err, error_msg);
}

/* align coords1 to coords2, no translation is done */
/* input: coords1(atom, dim), coords2(atom, dim), weight(atom) */
/* output: rotation(dim, dim), err */
CcpnStatus align_coordinates(int natoms, float **coords1, float **coords2,
                         float *weight, float **rotation, float *err, CcpnString error_msg)
{
    int i, j, k;
    float e, r, s, s0, s1, s2, w = 1;
    float work[COORD_NDIMS], eigenvalues[COORD_NDIMS], coords[COORD_NDIMS];
    static float RtR0[COORD_NDIMS], RtR1[COORD_NDIMS], RtR2[COORD_NDIMS];
    float *RtR[COORD_NDIMS] = { RtR0, RtR1, RtR2 };
    float **eigenvectors = RtR;
    float *e0 = eigenvectors[0], *e1 = eigenvectors[1], *e2 = eigenvectors[2];
    float v0[COORD_NDIMS], v1[COORD_NDIMS], v2[COORD_NDIMS];
    static float b0[COORD_NDIMS], b1[COORD_NDIMS], b2[COORD_NDIMS];
    float *b[] = { b0, b1, b2 };
    Bool find_eigenvectors = CCPN_TRUE;

    e = 0;
    for (i = 0; i < natoms; i++)
    {
        if (weight)
            w = weight[i];
        e += w * (inner_product(coords1[i], coords1[i], COORD_NDIMS) +
                  inner_product(coords2[i], coords2[i], COORD_NDIMS));
    }
    e *= 0.5;

    /*  calculate rotation */
    for (i = 0; i < COORD_NDIMS; i++)
    {
        for (j = 0; j < COORD_NDIMS; j++)
        {
            r = 0;
            for (k = 0; k < natoms; k++)
            {
                if (weight)
                    w = weight[k];
                r += w * coords2[k][i] * coords1[k][j];
            }
            rotation[i][j] = r;
        }
    }

    /* calculate RtR */
    for (i = 0; i < COORD_NDIMS; i++)
    {
        for (j = 0; j < COORD_NDIMS; j++)
        {
            r = 0;
            for (k = 0; k < COORD_NDIMS; k++)
                r += rotation[k][i] * rotation[k][j];
            RtR[i][j] = r;
        }
    }

/*
    print_matrix("rotation", rotation);
    print_matrix("RtR", RtR);
*/

    CHECK_STATUS(eigenvalues_eigenvectors(COORD_NDIMS, RtR, eigenvalues, work, find_eigenvectors, error_msg));

/*
    printf("evalues 1: %f %f %f\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
*/

    sort_eigenvalues(eigenvalues, eigenvectors);
/*
    printf("evalues 2: %f %f %f\n", eigenvalues[0], eigenvalues[1], eigenvalues[2]);
*/
    /* eigenvectors transpose */
    for (i = 0; i < COORD_NDIMS; i++)
    {
        for (j = i+1; j < COORD_NDIMS; j++)
            SWAP(eigenvectors[i][j], eigenvectors[j][i], float);
    }

    cross_product(e2, e0, e1);

    matrix_vector_multiply(v0, rotation, e0, COORD_NDIMS, COORD_NDIMS);
    matrix_vector_multiply(v1, rotation, e1, COORD_NDIMS, COORD_NDIMS);
    matrix_vector_multiply(v2, rotation, e2, COORD_NDIMS, COORD_NDIMS);

    s0 = (float) sqrt((double) eigenvalues[0]);
    s1 = (float) sqrt((double) eigenvalues[1]);
    s2 = (float) sqrt((double) eigenvalues[2]);

    s = 1 / s0;
    SCALE_VECTOR(b0, v0, s, COORD_NDIMS);
    s = 1 / s1;
    SCALE_VECTOR(b1, v1, s, COORD_NDIMS);
    cross_product(b2, b0, b1);

    for (i = 0; i < COORD_NDIMS; i++)
    {
        for (j = 0; j < COORD_NDIMS; j++)
            rotation[i][j] = b0[i]*e0[j] + b1[i]*e1[j] + b2[i]*e2[j];
    }

/*
    print_matrix("rotation", rotation);
*/

    if (inner_product(b2, v2, COORD_NDIMS) < 0)
        s = -1;
    else
        s = 1;

    *err = e - (s0 + s1 + s*s2);

    for (i = 0; i < natoms; i++)
    {   
        matrix_vector_multiply(coords, rotation, coords1[i], COORD_NDIMS, COORD_NDIMS);
        COPY_VECTOR(coords1[i], coords, COORD_NDIMS);
    }

    return CCPN_OK;
}

#define  RMSD_THRESHOLD  0.8

/* align ensemble of coords */
/* unlike align_coordinates, here weight cannot be null */
/* input: coords(ensemble, atom, dim), weight(atom) */
/* output: updated coords(ensemble, atom, dim), atom_rmsd(atom), ensemble_rmsd(ensemble), toterr */
CcpnStatus align_ensemble(int nensembles, int natoms, float ***coords,
			float *weight, float *atom_rmsd, float *ensemble_rmsd,
			float *toterr, CcpnString error_msg)
{
    int i, j, k, best_ensemble;
    float w, t, totw, err, delta, best_rmsd, mean[COORD_NDIMS];
    static float r0[COORD_NDIMS], r1[COORD_NDIMS], r2[COORD_NDIMS];
    float *rotation[] = {r0, r1, r2};

    /* first translate coords to center of mass */
    for (i = 0; i < nensembles; i++)
	translate_coordinates(natoms, coords[i], weight);

    /* fit against first member of ensemble */
    for (i = 1; i < nensembles; i++)
	CHECK_STATUS(align_coordinates(natoms, coords[i], coords[0], weight, rotation, &err, error_msg));

    /* calculate atom_rmsd and ensemble_rmsd */
    ZERO_VECTOR(ensemble_rmsd, nensembles);
    totw = 0;
    for (i = 0; i < natoms; i++)
    {
	w = weight[i];
	totw += w;

	ZERO_VECTOR(mean, COORD_NDIMS);
	for (j = 0; j < nensembles; j++)
	{
            for (k = 0; k < COORD_NDIMS; k++)
		mean[k] += coords[j][i][k];
	}

        for (k = 0; k < COORD_NDIMS; k++)
	    mean[k] /= nensembles;

	err = 0;
	for (j = 0; j < nensembles; j++)
	{
            for (k = 0; k < COORD_NDIMS; k++)
	    {
		delta = coords[j][i][k] - mean[k];
		delta *= delta;
		err += delta;
		ensemble_rmsd[j] += w*delta;
	    }
	}

	atom_rmsd[i] = sqrt(err / nensembles);
    }

    for (i = 0; i < nensembles; i++)
	ensemble_rmsd[i] = sqrt(ensemble_rmsd[i] / totw);

    /* determine ensemble closest to mean */
    best_ensemble = 0;
    best_rmsd = ensemble_rmsd[0];
    for (i = 1; i < nensembles; i++)
    {
	if (ensemble_rmsd[i] < best_rmsd)
	{
	    best_rmsd = ensemble_rmsd[i];
	    best_ensemble = i;
	}
    }

    /* adjust weight taking rmsd into account */
    for (i = 0; i < natoms; i++)
    {
	t = atom_rmsd[i] / RMSD_THRESHOLD;
	weight[i] *= exp(-t*t);
    }

    /* fit against best member of ensemble */
    *toterr = 0;
    for (i = 0; i < nensembles; i++)
    {
	if (i != best_ensemble)
	{
	    CHECK_STATUS(align_coordinates(natoms, coords[i], coords[best_ensemble], weight, rotation, &err, error_msg));
	    *toterr += err;
	}
    }

    /* calculate atom_rmsd and ensemble_rmsd again given RMSD weighted alignment*/
    ZERO_VECTOR(ensemble_rmsd, nensembles);
    totw = 0;
    for (i = 0; i < natoms; i++)
    {
	w = weight[i];
	totw += w;

	ZERO_VECTOR(mean, COORD_NDIMS);
	for (j = 0; j < nensembles; j++)
	{
            for (k = 0; k < COORD_NDIMS; k++)
		mean[k] += coords[j][i][k];
	}

        for (k = 0; k < COORD_NDIMS; k++)
	    mean[k] /= nensembles;

	err = 0;
	for (j = 0; j < nensembles; j++)
	{
            for (k = 0; k < COORD_NDIMS; k++)
	    {
		delta = coords[j][i][k] - mean[k];
		delta *= delta;
		err += delta;
		ensemble_rmsd[j] += w*delta;
	    }
	}

	atom_rmsd[i] = sqrt(err / nensembles);
    } 

    for (i = 0; i < nensembles; i++)
	ensemble_rmsd[i] = sqrt(ensemble_rmsd[i] / totw);

    return CCPN_OK;
}

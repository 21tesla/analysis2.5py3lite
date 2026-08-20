#include "../inc/myGEOMETRY.h"

int fuzzyEqual(const double d1, const double d2) {
    if ((d1 == 0.0) || (d2 == 0.0))
        return 1e-5 > fabs(d1 - d2);
    else
        return 1e-5 > fabs(1.0 - (d2 / d1));
}

/************************************************
 * Geometric Functions                          *
 ************************************************/

double AngleRestrict(const double theta) {
    /* synthetic extended precision constants */
    const double    P1 = 4 * 7.8539812564849853515625e-01;
    const double    P2 = 4 * 3.7748947079307981766760e-08;
    const double    P3 = 4 * 2.6951514290790594840552e-15;
    const double    TwoPi = 2 * (P1 + P2 + P3);

    const double    y = GSL_SIGN(theta) * 2 * floor(fabs(theta) / TwoPi);
    double          r = ((theta - y * P1) - y * P2) - y * P3;

    if (r > M_PI) {
        r = (((r - 2 * P1) - 2 * P2) - 2 * P3);
    } /* r-TwoPi */
    else if (r < -M_PI)
        r = (((r + 2 * P1) + 2 * P2) + 2 * P3); /* r+TwoPi */

    return r;
}

double NormNormalize(gsl_vector * vector) {
    double          norm = Norm(vector);

    gsl_vector_scale(vector, 1. / norm);
    return norm;
}

void CrossProduct(const gsl_vector * first, const gsl_vector * second, gsl_vector * result) {
    gsl_vector_set(result, 0, Y(first) * Z(second) - Y(second) * Z(first));
    gsl_vector_set(result, 1, Z(first) * X(second) - Z(second) * X(first));
    gsl_vector_set(result, 2, X(first) * Y(second) - X(second) * Y(first));
}

double Angle(const gsl_vector * first, const gsl_vector * second) {
    double          n1 = Norm(first), n2 = Norm(second);
    double          c;

    DotProduct(first, second, &c);

    if (n1 * n2)
        c /= n1 * n2;
    else
        return PI;

    c = GSL_MAX_DBL(c, -1.);
    c = GSL_MIN_DBL(c, 1.);

    return acos(c);
}

void Points2Vector(const gsl_vector * p1, const gsl_vector * p2, gsl_vector * v) {
    gsl_vector_memcpy(v, p2);
    gsl_vector_sub(v, p1);
}

double Points2UVector(const gsl_vector * p1, const gsl_vector * p2, gsl_vector * v) {
    gsl_vector_memcpy(v, p2);
    gsl_vector_sub(v, p1);
    return NormNormalize(v);
}

double Dist(const gsl_vector * p1, const gsl_vector * p2) {
    double          dist;
    gsl_vector     *v = gsl_vector_alloc(p1->size);

    Points2Vector(p1, p2, v);
    dist = Norm(v);

    gsl_vector_free(v);
    return dist;
}

double CalcAngle(const gsl_vector * p1, const gsl_vector * p2, const gsl_vector * p3) {
    double          angle;
    gsl_vector     *v1 = gsl_vector_alloc(3);
    gsl_vector     *v2 = gsl_vector_alloc(3);

    Points2Vector(p2, p1, v1);
    Points2Vector(p2, p3, v2);

    angle = Angle(v1, v2);

    gsl_vector_free(v1);
    gsl_vector_free(v2);

    return angle;
}

double CalcDihedral(const gsl_vector * p1, const gsl_vector * p2, const gsl_vector * p3, const gsl_vector * p4) {
    double          angle;

    gsl_vector     *ab = gsl_vector_alloc(3);
    gsl_vector     *cb = gsl_vector_alloc(3);
    gsl_vector     *db = gsl_vector_alloc(3);

    gsl_vector     *u = gsl_vector_alloc(3);
    gsl_vector     *v = gsl_vector_alloc(3);
    gsl_vector     *w = gsl_vector_alloc(3);

    Points2Vector(p2, p1, ab);
    Points2Vector(p2, p3, cb);
    Points2Vector(p3, p4, db);

    CrossProduct(ab, cb, u);
    CrossProduct(db, cb, v);
    CrossProduct(u, v, w);

    angle = Angle(u, v);

    if (Angle(cb, w) > 0.001)
        angle = -angle;

    gsl_vector_free(ab);
    gsl_vector_free(cb);
    gsl_vector_free(db);

    gsl_vector_free(u);
    gsl_vector_free(v);
    gsl_vector_free(w);

    return angle;
}

void CalcPhiPsiTet(const gsl_vector * c_n__1, const gsl_vector * n_ca_1, const gsl_vector * ca_c_2, const gsl_vector * c_n__2, double *phi, double *psi, double *tet) {
    gsl_vector     *u_phi = gsl_vector_alloc(3);
    gsl_vector     *v_phi = gsl_vector_alloc(3);
    gsl_vector     *w_phi = gsl_vector_alloc(3);

    gsl_vector     *u_psi = gsl_vector_alloc(3);
    gsl_vector     *v_psi = gsl_vector_alloc(3);
    gsl_vector     *w_psi = gsl_vector_alloc(3);

    CrossProduct(n_ca_1, c_n__1, u_phi);
    CrossProduct(ca_c_2, n_ca_1, v_phi);
    CrossProduct(u_phi, v_phi, w_phi);

    gsl_vector_memcpy(u_psi, v_phi);
    CrossProduct(c_n__2, ca_c_2, v_psi);
    CrossProduct(u_psi, v_psi, w_psi);

    (*phi) = Angle(u_phi, v_phi);
    if (Angle(n_ca_1, w_phi) > 0.001)
        (*phi) = -(*phi);

    (*psi) = Angle(u_psi, v_psi);
    if (Angle(ca_c_2, w_psi) > 0.001)
        (*psi) = -(*psi);

    (*tet) = PI - Angle(n_ca_1, ca_c_2);

    gsl_vector_free(u_phi);
    gsl_vector_free(v_phi);
    gsl_vector_free(w_phi);

    gsl_vector_free(u_psi);
    gsl_vector_free(v_psi);
    gsl_vector_free(w_psi);
}

void RotMat(const gsl_vector * u, double theta, gsl_matrix * rotMat) {
    double          c = cos(theta);
    double          s = sin(theta);
    double          t = 1 - c;
    double          norm = Norm(u);
    double          x = X(u) / norm;
    double          y = Y(u) / norm;
    double          z = Z(u) / norm;

    // first Row
    gsl_matrix_set(rotMat, 0, 0, t * x * x + c);
    gsl_matrix_set(rotMat, 0, 1, t * x * y - s * z);
    gsl_matrix_set(rotMat, 0, 2, t * x * z + s * y);

    // Second Row
    gsl_matrix_set(rotMat, 1, 0, t * x * y + s * z);
    gsl_matrix_set(rotMat, 1, 1, t * y * y + c);
    gsl_matrix_set(rotMat, 1, 2, t * y * z - s * x);

    // Third Row
    gsl_matrix_set(rotMat, 2, 0, t * x * z - s * y);
    gsl_matrix_set(rotMat, 2, 1, t * y * z + s * x);
    gsl_matrix_set(rotMat, 2, 2, t * z * z + c);
}

void AxisAngle2Mat(const gsl_vector * u, gsl_matrix * rotMat) {
    double          theta = Norm(u);

    if (theta == 0.) {
        gsl_matrix_set_identity(rotMat);
    } else {
        double          x = X(u) / theta;
        double          y = Y(u) / theta;
        double          z = Z(u) / theta;
        double          c = cos(theta);
        double          s = sin(theta);
        double          t = 1 - c;
        double          tmp1, tmp2;

        gsl_matrix_set(rotMat, 0, 0, t * x * x + c);
        gsl_matrix_set(rotMat, 1, 1, t * y * y + c);
        gsl_matrix_set(rotMat, 2, 2, t * z * z + c);

        tmp1 = t * x * y;
        tmp2 = s * z;
        gsl_matrix_set(rotMat, 0, 1, tmp1 - tmp2);
        gsl_matrix_set(rotMat, 1, 0, tmp1 + tmp2);

        tmp1 = t * x * z;
        tmp2 = s * y;
        gsl_matrix_set(rotMat, 0, 2, tmp1 + tmp2);
        gsl_matrix_set(rotMat, 2, 0, tmp1 - tmp2);

        tmp1 = t * y * z;
        tmp2 = s * x;
        gsl_matrix_set(rotMat, 1, 2, tmp1 - tmp2);
        gsl_matrix_set(rotMat, 2, 1, tmp1 + tmp2);
    }
}

void Mat2AxisAngle(const gsl_matrix * rotMat, gsl_vector * u) {
    double          angle;
    double          cosTheta = (gsl_matrix_get(rotMat, 0, 0) + gsl_matrix_get(rotMat, 1, 1) + gsl_matrix_get(rotMat, 2, 2) - 1.0) / 2.0;

    if (fuzzyEqual(cosTheta, 1)) {
        gsl_vector_set_zero(u);
        return;
    }
    if (!fuzzyEqual(cosTheta, -1.0)) {
        double          sinTheta = sqrt(1.0 - cosTheta * cosTheta);

        angle = acos(cosTheta);
        gsl_vector_set(u, 0, (gsl_matrix_get(rotMat, 2, 1) - gsl_matrix_get(rotMat, 1, 2)) / 2.0 / sinTheta);
        gsl_vector_set(u, 1, (gsl_matrix_get(rotMat, 0, 2) - gsl_matrix_get(rotMat, 2, 0)) / 2.0 / sinTheta);
        gsl_vector_set(u, 2, (gsl_matrix_get(rotMat, 1, 0) - gsl_matrix_get(rotMat, 0, 1)) / 2.0 / sinTheta);
        gsl_vector_scale(u, angle);
        return;
    }
    angle = PI;
    if (gsl_matrix_get(rotMat, 0, 0) >= gsl_matrix_get(rotMat, 1, 1)) {
        if (gsl_matrix_get(rotMat, 0, 0) >= gsl_matrix_get(rotMat, 2, 2)) {
            // rot00 is maximal diagonal term
            double          halfInverse = 0.0;

            gsl_vector_set(u, 0, sqrt(gsl_matrix_get(rotMat, 0, 0) - gsl_matrix_get(rotMat, 1, 1) - gsl_matrix_get(rotMat, 2, 2) + 1.0) / 2.0);
            halfInverse = 1.0 / (2.0 * gsl_vector_get(u, 0));
            gsl_vector_set(u, 1, halfInverse * gsl_matrix_get(rotMat, 0, 1));
            gsl_vector_set(u, 2, halfInverse * gsl_matrix_get(rotMat, 0, 2));
            gsl_vector_scale(u, angle);
        } else {
            // rot22 is maximal diagonal term
            double          halfInverse = 0.0;

            gsl_vector_set(u, 2, sqrt(gsl_matrix_get(rotMat, 2, 2) - gsl_matrix_get(rotMat, 0, 0) - gsl_matrix_get(rotMat, 1, 1) + 1.0) / 2.0);
            halfInverse = 1.0 / (2.0 * gsl_vector_get(u, 2));
            gsl_vector_set(u, 0, halfInverse * gsl_matrix_get(rotMat, 0, 2));
            gsl_vector_set(u, 1, halfInverse * gsl_matrix_get(rotMat, 1, 2));
            gsl_vector_scale(u, angle);
        }
    } else {
        if (gsl_matrix_get(rotMat, 1, 1) >= gsl_matrix_get(rotMat, 2, 2)) {
            // rot11 is maximal diagonal term
            double          halfInverse = 0.0;

            gsl_vector_set(u, 1, sqrt(gsl_matrix_get(rotMat, 1, 1) - gsl_matrix_get(rotMat, 0, 0) - gsl_matrix_get(rotMat, 2, 2) + 1.0) / 2.0);
            halfInverse = 1.0 / (2.0 * gsl_vector_get(u, 1));
            gsl_vector_set(u, 0, halfInverse * gsl_matrix_get(rotMat, 0, 1));
            gsl_vector_set(u, 2, halfInverse * gsl_matrix_get(rotMat, 1, 2));
            gsl_vector_scale(u, angle);
        } else {
            // rot22 is maximal diagonal term
            double          halfInverse = 0.0;

            gsl_vector_set(u, 2, sqrt(gsl_matrix_get(rotMat, 2, 2) - gsl_matrix_get(rotMat, 0, 0) - gsl_matrix_get(rotMat, 1, 1) + 1.0) / 2.0);
            halfInverse = 1.0 / (2.0 * gsl_vector_get(u, 2));
            gsl_vector_set(u, 0, halfInverse * gsl_matrix_get(rotMat, 0, 2));
            gsl_vector_set(u, 1, halfInverse * gsl_matrix_get(rotMat, 1, 2));
            gsl_vector_scale(u, angle);
        }
    }
}

/*======================================================================*
 *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991               *
 *                                                                      *
 *  This routine maps three values (x[0], x[1], x[2]) in the range      *
 *  [0,1] into a 3x3 rotation matrix M.  Uniformly distributed random   *
 *  variables x0, x1, and x2 create uniformly distributed random        *
 *  rotation matrices.  To create small uniformly distributed           *
 *  "perturbations", supply samples in the following ranges             *
 *                                                                      *
 *      x[0] in [ 0, d ]                                                *
 *      x[1] in [ 0, 1 ]                                                *
 *      x[2] in [ 0, d ]                                                *
 *                                                                      *
 * where 0 < d < 1 controls the size of the perturbation.  Any of the   *
 * random variables may be stratified (or "jittered") for a slightly    *
 * more even distribution.                                              *
 *                                                                      *
 *======================================================================*/

void RandMat(const gsl_vector * u, gsl_matrix * rotMat) {

    double          theta = X(u) * 2. * PI; /* Rotation about the pole (Z).  */
    double          phi = Y(u) * 2. * PI;   /* For direction of pole deflection. */
    double          z = Z(u) * 2.;  /* For magnitude of pole deflection. */

    /* Compute a vector V used for distributing points over the sphere */
    /* via the reflection I - V Transpose(V).  This formulation of V */
    /* will guarantee that if x[1] and x[2] are uniformly distributed, */
    /* the reflected points will be uniform on the sphere.  Note that V */
    /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */

    double          r = sqrt(z);
    double          Vx = sin(phi) * r;
    double          Vy = cos(phi) * r;
    double          Vz = sqrt(2.0 - z);

    /* Compute the row vector S = Transpose(V) * R, where R is a simple */
    /* rotation by theta about the z-axis.  No need to compute Sz since */
    /* it's just Vz.  */

    double          st = sin(theta);
    double          ct = cos(theta);
    double          Sx = Vx * ct - Vy * st;
    double          Sy = Vx * st + Vy * ct;

    /* Construct the rotation matrix ( V Transpose(V) - I ) R, which */
    /* is equivalent to V S - R.  */

    // first Row
    gsl_matrix_set(rotMat, 0, 0, Vx * Sx - ct);
    gsl_matrix_set(rotMat, 0, 1, Vx * Sy - st);
    gsl_matrix_set(rotMat, 0, 2, Vx * Vz);

    // Second Row
    gsl_matrix_set(rotMat, 1, 0, Vy * Sx + st);
    gsl_matrix_set(rotMat, 1, 1, Vy * Sy - ct);
    gsl_matrix_set(rotMat, 1, 2, Vy * Vz);

    // Third Row
    gsl_matrix_set(rotMat, 2, 0, Vz * Sx);
    gsl_matrix_set(rotMat, 2, 1, Vz * Sy);
    gsl_matrix_set(rotMat, 2, 2, 1.0 - z);  /* This equals Vz * Vz - 1.0 */
}

void RotationTranslation(gsl_matrix * matrix, gsl_vector * transl, gsl_vector * vector, gsl_vector * result) {
    gsl_vector_memcpy(result, transl);
    gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, vector, 1.0, result);
}

void ScalingTranslation(double scale, gsl_vector * transl, gsl_vector * vector, gsl_vector * result) {
    gsl_vector_memcpy(result, vector);
    gsl_vector_scale(result, scale);
    gsl_vector_add(result, transl);
}

void Cart2Spher(const gsl_vector * vector, double *rho, double *theta, double *phi) {
    double          vx = gsl_vector_get(vector, 0);
    double          vy = gsl_vector_get(vector, 1);
    double          vz = gsl_vector_get(vector, 2);

    (*rho) = gsl_blas_dnrm2(vector);

    if (vx && (*rho)) {
        *theta = acos(vz / (*rho));
        *phi = atan(vy / vx);
        if (vx < 0)
            (*phi) += PI;
    } else {
        *theta = acos(vz);
        *phi = (double) GSL_SIGN(vy) * M_PI_2;
    }
}

void RotEuler(const gsl_vector * vector, gsl_vector * result, double alpha, double beta, double gamma) {
    double          ca = cos(alpha), sa = sin(alpha);
    double          cb = cos(beta), sb = sin(beta);
    double          cg = cos(gamma), sg = sin(gamma);

    double          ROT[] = {
        ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb,
        sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg, sa * sb,
        -sb * cg, sb * sg, cb
    };

    gsl_matrix_view rot = gsl_matrix_view_array(ROT, 3, 3);

    Rotation(&rot.matrix, vector, result);
}

gsl_matrix     *MatEuler(double alpha, double beta, double gamma) {
    gsl_matrix     *euler = gsl_matrix_calloc(3, 3);
    double          ca = cos(alpha), sa = sin(alpha);
    double          cb = cos(beta), sb = sin(beta);
    double          cg = cos(gamma), sg = sin(gamma);

    double          ROT[] = {
        ca * cg - sa * cb * sg, -sa * cg - ca * cb * sg, sb * sg,
        ca * sg + sa * cb * cg, -sa * sg + ca * cb * cg, -sb * cg,
        sa * sb, ca * sb, cb
    };

    gsl_matrix_view rot = gsl_matrix_view_array(ROT, 3, 3);

    gsl_matrix_memcpy(euler, &rot.matrix);

    return euler;
}

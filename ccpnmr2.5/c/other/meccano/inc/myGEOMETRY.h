#ifndef _MYGEO_H_
#define _MYGEO_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#define PI M_PI

#define RAD2DEG(angle) (180./PI*(angle))
#define DEG2RAD(angle) (PI/180.*(angle))

#define X(vector) gsl_vector_get(vector,0)
#define Y(vector) gsl_vector_get(vector,1)
#define Z(vector) gsl_vector_get(vector,2)

#define Norm(vector) gsl_blas_dnrm2(vector)
#define Normalize(vector) gsl_vector_scale(vector, 1./gsl_blas_dnrm2(vector))
#define DotProduct(first, second, result) gsl_blas_ddot(first, second, result)

#define Rotation(matrix, vector, result) gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, vector, 0.0, result)
#define MatrixMultiplication(matrix1, matrix2, result) gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix1, matrix2, 0.0, result)

extern double   AngleRestrict(const double theta);

extern double   NormNormalize(gsl_vector * vector);

extern void     CrossProduct(const gsl_vector * first, const gsl_vector * second, gsl_vector * result);

extern double   Angle(const gsl_vector * first, const gsl_vector * second);

extern void     Points2Vector(const gsl_vector * p1, const gsl_vector * p2, gsl_vector * v);

extern double   Points2UVector(const gsl_vector * p1, const gsl_vector * p2, gsl_vector * v);

extern double   Dist(const gsl_vector * p1, const gsl_vector * p2);

extern double   CalcAngle(const gsl_vector * p1, const gsl_vector * p2, const gsl_vector * p3);

extern double   CalcDihedral(const gsl_vector * p1, const gsl_vector * p2, const gsl_vector * p3, const gsl_vector * p4);

extern void     CalcPhiPsiTet(const gsl_vector * c_n__1, const gsl_vector * n_ca_1, const gsl_vector * ca_c_2, const gsl_vector * c_n__2, double *phi, double *psi, double *tet);

extern void     RotMat(const gsl_vector * u, double theta, gsl_matrix * rotMat);

extern void     AxisAngle2Mat(const gsl_vector * u, gsl_matrix * rotMat);
extern void     Mat2AxisAngle(const gsl_matrix * rotMat, gsl_vector * u);

extern void     RandMat(const gsl_vector * u, gsl_matrix * rotmat);

extern void     RotationTranslation(gsl_matrix * matrix, gsl_vector * transl, gsl_vector * vector, gsl_vector * result);

extern void     ScalingTranslation(double scale, gsl_vector * transl, gsl_vector * vector, gsl_vector * result);

extern void     Cart2Spher(const gsl_vector * vector, double *rho, double *theta, double *phi);

extern void     RotEuler(const gsl_vector * vector, gsl_vector * result, double alpha, double beta, double gamma);

extern gsl_matrix *MatEuler(double alpha, double beta, double gamma);
#endif

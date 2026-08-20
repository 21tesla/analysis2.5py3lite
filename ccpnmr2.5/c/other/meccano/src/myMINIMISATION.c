#include "../inc/myMINIMISATION.h"

extern gsl_rng *r;
extern int      printFlag;


/*
 * PosA <-- PosB
 */
void PositionCpy(position_t * posA, position_t * posB) {
    posA->f = posB->f;
    gsl_vector_memcpy(posA->x, posB->x);
}


/*
 * PSO Minimization
 */
double PSOMinRotFwd(double (*Func) (const gsl_vector *, void *), gsl_vector * x, void *params) {
    int             dimNb;      // Search space dimension
    int             d;          // Current dimensionn
    int             swarmSize;  // Swarm size
    int             s;          // Rank of the current particle
    int             ss;
    int             g;          // Rank of the best informant
    int             best;       // Best of the best position (rank in the swarm)
    position_t     *position;   // Positions
    position_t     *bestPosition;   // Best positions found by each particle
    position_t      avPos;      // Barycentre
    gsl_vector     *xMin, *xMax;    // Intervals defining the search space
    gsl_matrix_int *links;      // Information links
    int             initLinks;  // Flag to (re)init or not the information links
    int             maxNbPartInf;   // Max number of particles informed by a given one

    double          w;          // first cidence coefficient

    int             nbEval;     // Total number of evaluations
    int             nbEvalMax;  // Max number of evaluations

    int             i;
    double          chi2Best, chi2BestPrev;
    double          mean, stdev;

    int             ppi = ((param_t *) params)->ppi;
    int             pPlaneNbInMin = ((param_t *) params)->pPlaneNbInMin;
    int             windowStart = ((param_t *) params)->windowStart;
    polyPeptide_t  *polyP = ((param_t *) params)->polyP;
    ramaDB_t       *ramaDB = ((param_t *) params)->ramaDB;

    double          spaceSize;

    int             k = 0;
    double          diversity = 1e32;
    double phi, psi;

    int             j;

    // ---------------- PARAMETERS
    dimNb = x->size;            // Search space dimension
    swarmSize = (10 + (int) (2. * sqrt(dimNb)));
    maxNbPartInf = 3;
    w = 0.74 /* 1. / (2. * log(2.)) */ ;

    nbEvalMax = (1000 * swarmSize);    // Max number of evaluations for each

    // printf("\n Swarm size %i", S);
    // printf("\n coefficients %f %f \n", w, c);
    // -----------------------------------------------------

    // INITIALISATION

    position = (position_t *) malloc(swarmSize * sizeof(position_t));   // Positions
    bestPosition = (position_t *) malloc(swarmSize * sizeof(position_t));   // Best positions found by each particle
    for (s = 0; s < swarmSize; s++) {
        position[s].x = gsl_vector_alloc(dimNb);
        bestPosition[s].x = gsl_vector_alloc(dimNb);
    }
    avPos.x = gsl_vector_alloc(dimNb);  // Barycentre
    xMin = gsl_vector_alloc(dimNb);
    xMax = gsl_vector_alloc(dimNb); // Intervals defining the search space
    links = gsl_matrix_int_alloc(swarmSize, swarmSize); // Information links

    // D-cube data_t
    spaceSize = 0.;
    for (j = 0; j < pPlaneNbInMin; j++) {
        if (windowStart != 0 && !j) {
            gsl_vector_view xMinj = gsl_vector_subvector(xMin, 3 * j, 3);
            gsl_vector_view xMaxj = gsl_vector_subvector(xMax, 3 * j, 3);

            SET_ALL(&xMinj.vector, -PI);
            SET_ALL(&xMaxj.vector, PI);
            spaceSize += 3.0 * gsl_pow_2(2.0 * PI);
        } else {
            gsl_vector_view xMinj = gsl_vector_subvector(xMin, 3 * j, 3);
            gsl_vector_view xMaxj = gsl_vector_subvector(xMax, 3 * j, 3);

            SET(&xMinj.vector, 0, -PI);
            SET(&xMaxj.vector, 0, PI);
            spaceSize += gsl_pow_2(2.0 * PI);
            SET(&xMinj.vector, 1, -PI);
            SET(&xMaxj.vector, 1, PI);
            spaceSize += gsl_pow_2(2.0 * PI);
            SET(&xMinj.vector, 2, DEG2RAD(95.0));
            SET(&xMaxj.vector, 2, DEG2RAD(125.0));
            spaceSize += gsl_pow_2(DEG2RAD(30.0));
        }
    }

    spaceSize = sqrt(spaceSize);

    // Initialisation of information variables
    for (s = 0; s < swarmSize; s++) // Positions and velocities
    {
        for (j = 0; j < pPlaneNbInMin; j++) {
            if (windowStart != 0 && !j) {
                gsl_vector     *u = gsl_vector_alloc(3);
                gsl_matrix     *rotMat = gsl_matrix_alloc(3, 3);
                gsl_vector_view xView = gsl_vector_subvector(position[s].x, 3 * j, 3);

                SET(u, 0, ALEA(0.0, 1.0));
                SET(u, 1, ALEA(0.0, 1.0));
                SET(u, 2, ALEA(0.0, 1.0));
                RandMat(u, rotMat);
                Mat2AxisAngle(rotMat, &xView.vector);

                gsl_vector_free(u);
                gsl_matrix_free(rotMat);
            } else {
                AleaRamaDB(ramaDB, polyP, ppi + j - 1, &phi, &psi);
                SET(position[s].x, 0 + 3 * j, phi);
                SET(position[s].x, 1 + 3 * j, psi);
                SET(position[s].x, 2 + 3 * j, gsl_ran_gaussian(r, polyP->dTet0) + polyP->tet0);
            }
        }
    }

    // first evaluations
    nbEval = 0;
    for (s = 0; s < swarmSize; s++) {
        nbEval++;
        position[s].f = (*Func) (position[s].x, params);
    }

    for (s = 0; s < swarmSize; s++) {
        PositionCpy(&bestPosition[s], &position[s]);    // Best position = current one
    }

    // Find the best
    best = 0;
    for (s = 0; s < swarmSize; s++) {
        if (bestPosition[s].f < bestPosition[best].f) {
            best = s;
        }
    }
    chi2Best = bestPosition[best].f;       // Current min chi2Best
    chi2BestPrev = chi2Best;    // Previous min chi2Best

    initLinks = 1;              // So that information links will be initialized

    // ---------------------------------------------- ITERATIONS

    while (diversity > 1e-5 && nbEval < nbEvalMax) {

        if (initLinks == 1) {   // Who informs who, at random
            for (s = 0; s < swarmSize; s++) {
                for (ss = 0; ss < swarmSize; ss++) {
                    IMSET(links, ss, s, 0); // Init to "no link"
                }
                IMSET(links, s, s, 1);  // Each particle informs itself
            }

            for (ss = 0; ss < swarmSize; ss++)  // Other links
            {
                for (i = 0; i < maxNbPartInf; i++) {
                    s = IALEA(0, swarmSize - 1);
                    IMSET(links, ss, s, 1);
                }
            }
        }
        // The swarm MOVES
        for (s = 0; s < swarmSize; s++) // For each particle ...
        {
            // .. find the best informant
            g = s;
            for (ss = 0; ss < swarmSize; ss++) {
                if (IMGET(links, ss, s) == 1 && bestPosition[ss].f < bestPosition[g].f)
                    g = ss;
            }

            // ... compute the new velocity, and move

            for (d = 0; d < dimNb; d++) {
                mean = (GET(bestPosition[s].x, d) + GET(bestPosition[g].x, d)) / 2.;
                stdev = w * fabs(GET(bestPosition[s].x, d) - GET(bestPosition[g].x, d));
                SET(position[s].x, d, mean + gsl_ran_gaussian(r, stdev));
            }

            // ... interval cinement (keep in the box)
            for (d = 0; d < dimNb; d++) {
                if (GET(position[s].x, d) < GET(xMin, d)) {
                    SET(position[s].x, d, GET(xMin, d));
                }
                if (GET(position[s].x, d) > GET(xMax, d)) {
                    SET(position[s].x, d, GET(xMax, d));
                }
            }
            // ... evaluate the new position
            nbEval++;
            position[s].f = (*Func) (position[s].x, params);

            // ... update the best previous position
            if (position[s].f < bestPosition[s].f) {
                PositionCpy(&bestPosition[s], &position[s]);
                // ... update the best of the bests
                if (bestPosition[s].f < bestPosition[best].f)
                    best = s;
            }

            // Check if finished
            // If no improvement, information links will beta reinitialized
            chi2Best = bestPosition[best].f;
            if (chi2Best >= chi2BestPrev)
                initLinks = 1;
            else
                initLinks = 0;
            chi2BestPrev = chi2Best;

            if (k % 5 == 0) {
                diversity = 0.;
                gsl_vector_set_zero(avPos.x);

                for (s = 0; s < swarmSize; s++)
                    gsl_vector_add(avPos.x, position[s].x);
                gsl_vector_scale(avPos.x, 1.0 / (double) swarmSize);

                for (s = 0; s < swarmSize; s++)
                    diversity += Dist(position[s].x, avPos.x);
                diversity /= (swarmSize * spaceSize);
            }
            k++;
        }
    }

    if (printFlag)
        printf("PSO: %7d %.1e %lf\n", nbEval, diversity, bestPosition[best].f);

    gsl_vector_memcpy(x, bestPosition[best].x);

    gsl_matrix_int_free(links);
    gsl_vector_free(xMin);
    gsl_vector_free(xMax);
    gsl_vector_free(avPos.x);  // Barycentre
    for (s = 0; s < swarmSize; s++) {
        gsl_vector_free(position[s].x);
        gsl_vector_free(bestPosition[s].x);
    }
    free(bestPosition);
    free(position);

    return bestPosition[best].f;
}

// double PSOMinRotRev(double (*Func) (const gsl_vector *, void *), gsl_vector * x, void *params, ramaDB_t * ramaDB, polyPeptide_t * polyP, int ppi) {
//     int             best;       // Best of the best position (rank in the swarm)
//     int             D;          // Search space dimension
//     int             LINKS[S_max][S_max];    // Information links
//     position_t      P[S_max];   // Best positions found by each particle
//     int             S;          // Swarm size
//     position_t      X[S_max];   // Positions
//     double          xmin[D_max], xmax[D_max];   // Intervals defining the search space
//     double          w;          // first confidence coefficient
//     double          c;          // Second confidence coefficient
//     int             d;          // Current dimensionn
//     int             nb_eval;    // Total number of evaluations
//     int             eval_max;   // Max number of evaluations
//     int             g;          // Rank of the best informant
//     int             init_links; // Flag to (re)init or not the information links
//     int             i;
//     int             K;          // Max number of particles informed by a given one
//     int             s;          // Rank of the current particle
//     int             m;
//     double          chi2Best, chi2BestPrev;
//     double          mean, stdev;
// 
//     gsl_vector     *var;        // For GSL functions
//     double          L;
//     int             k = 0;
//     position_t      Av;         // Barycentre
//     double          diversity = 1e32;
// 
//     D = x->size;                // Search space dimension
//     var = gsl_vector_alloc(D);
//     double          phi, psi;
//     int             ii;
// 
//     // D-cube data_t
//     L = 0.;
//     for (d = 0; d < D; d++) {
//         xmin[d] = 0.;
//         xmax[d] = 1.;
//         L += 1.;
//     }
//     L = sqrt(L);
// 
//     eval_max = (10000 * (10 + (int) (2. * sqrt(D)))); // Max number of
//     // evaluations for each
//     // run
// 
//     // ---------------- PARAMETERS
//     S = (10 + (int) (2. * sqrt(D)));
//     if (S > S_max)
//         S = S_max;
//     K = 3;
//     w = 0.74 /* 1. / (2. * log(2.)) */ ;
//     c = 0.5 + log(2.);
//     // printf("\n Swarm size %i", S);
//     // printf("\n coefficients %f %f \n", w, c);
//     // -----------------------------------------------------
// 
//     // INITIALISATION
//     // Initialisation of information variables
// 
//     for (s = 0; s < S; s++)     // Positions and velocities
//     {
//         for (ii = 0; ii < D / 3; ii++) {
//             if (polyP->firstPepPlaNb && !ii) {
//                 X[s].x[3 * ii + 0] = ALEA(0., 1.);
//                 X[s].x[3 * ii + 1] = ALEA(0., 1.);
//                 X[s].x[3 * ii + 2] = ALEA(0., 1.);
//             } else {
//                 AleaRamaDB(ramaDB, polyP, ppi - ii, &phi, &psi);
//                 X[s].x[3 * ii + 0] = phi;
//                 X[s].x[3 * ii + 1] = psi;
//                 X[s].x[3 * ii + 2] = gsl_ran_gaussian(r, polyP->dTet0) + polyP->tet0;
// //                 printf("%8.3lf %8.3lf %8.3lf\n",X[s].x[3*ii+0], X[s].x[3*ii+1], X[s].x[3*ii+2]);
//             }
//         }
//     }
// 
//     // first evaluations
//     nb_eval = 0;
// 
//     for (s = 0; s < S; s++) {
//         nb_eval++;
//         for (i = 0; i < D; i++)
//             SET(var, i, X[s].x[i]);
//         X[s].f = (*Func) (var, params);
//         P[s] = X[s];            // Best position = current one
//     }
// 
//     // Find the best
//     best = 0;
//     for (s = 1; s < S; s++)
//         if (P[s].f < P[best].f)
//             best = s;
//     chi2Best = P[best].f;       // Current min chi2Best
//     chi2BestPrev = chi2Best;    // Previous min chi2Best
// 
//     init_links = 1;             // So that information links will be
//     // initialized
// 
//     // ---------------------------------------------- ITERATIONS
// 
//     // if (nb_eval < eval_max)
//     // goto loop;
//     while (diversity > 1e-5 && nb_eval < eval_max) {
// 
//         if (init_links == 1) {  // Who informs who, at random
//             for (s = 0; s < S; s++) {
//                 for (m = 0; m < S; m++)
//                     LINKS[m][s] = 0;    // Init to "no link"
//                 LINKS[s][s] = 1;    // Each particle informs itself
//             }
// 
//             for (m = 0; m < S; m++) // Other links
//             {
//                 for (i = 0; i < K; i++) {
//                     s = IALEA(0, S - 1);
//                     LINKS[m][s] = 1;
//                 }
//             }
//         }
//         // The swarm MOVES
//         for (s = 0; s < S; s++) // For each particle ...
//         {
//             // .. find the best informant
//             g = s;
//             for (m = 0; m < S; m++) {
//                 if (LINKS[m][s] == 1 && P[m].f < P[g].f)
//                     g = m;
//             }
// 
//             // ... compute the new velocity, and move
// 
//             for (d = 0; d < D; d++) {
//                 mean = (P[s].x[d] + P[g].x[d]) / 2.;
//                 stdev = w * fabs(P[s].x[d] - P[g].x[d]);
//                 X[s].x[d] = mean + gsl_ran_gaussian(r, stdev);
//             }
// 
//             // ... interval confinement (keep in the box)
//             for (d = 0; d < D; d++) {
//                 if (X[s].x[d] < xmin[d]) {
//                     X[s].x[d] = xmin[d];
//                 }
//                 if (X[s].x[d] > xmax[d]) {
//                     X[s].x[d] = xmax[d];
//                 }
//             }
// 
//             // ... evaluate the new position
//             nb_eval++;
//             for (i = 0; i < D; i++)
//                 SET(var, i, X[s].x[i]);
//             X[s].f = (*Func) (var, params);
// 
//             // ... update the best previous position
//             if (X[s].f < P[s].f) {
//                 P[s] = X[s];
//                 // ... update the best of the bests
//                 if (P[s].f < P[best].f)
//                     best = s;
//             }
//         }
// 
//         if (k % 5 == 0) {
//             for (d = 0; d < D; d++) {
//                 Av.x[d] = 0.;
//                 for (s = 0; s < S; s++)
//                     Av.x[d] += X[s].x[d];
//                 Av.x[d] /= S;
//             }
//             diversity = 0.;
//             for (s = 0; s < S; s++) {
//                 double          norm2 = 0.;
// 
//                 for (d = 0; d < D; d++)
//                     norm2 += gsl_pow_2(X[s].x[d] - Av.x[d]);
//                 diversity += sqrt(norm2);
//             }
//             diversity /= (S * L);
//         }
//         k++;
// 
//         // Check if finished
//         // If no improvement, information links will beta reinitialized
//         chi2Best = P[best].f;
//         if (chi2Best >= chi2BestPrev)
//             init_links = 1;
//         else
//             init_links = 0;
//         chi2BestPrev = chi2Best;
//     }
// 
//     if (printFlag)
//         printf("PSO: %7d %.1e %lf\n", nb_eval, diversity, P[best].f);
// 
//     for (i = 0; i < D; i++)
//         SET(x, i, P[best].x[i]);
// 
//     return P[best].f;
// }

double ConjGradMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params, int iterNb) {
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf f;
    int             status;
    int             iter = 0;
    int             p = x->size;

    f.f = Func_f;
    f.df = Func_df;
    f.fdf = Func_fdf;
    f.n = p;
    f.params = params;

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc(T, p);

    gsl_multimin_fdfminimizer_set(s, &f, x, 1e-5, 1e-3);

    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status)
            break;

        if (printFlag && iter % 50 == 0)
            printf("ConjGradMin: %3d    %8.3lf     %8.3lf\n", iter, s->f, gsl_blas_dnrm2(s->gradient));

        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
    } while (status == GSL_CONTINUE && iter < iterNb);

    if (printFlag)
        printf("ConjGradMin: %3d    %8.3lf     %8.3lf\n", iter, s->f, gsl_blas_dnrm2(s->gradient));

    gsl_vector_memcpy(x, s->x);

    gsl_multimin_fdfminimizer_free(s);

    return s->f;
}

double SteepDescMin(double (*Func_f) (const gsl_vector *, void *), void (*Func_df) (const gsl_vector *, void *, gsl_vector *), void (*Func_fdf) (const gsl_vector *, void *, double *, gsl_vector *), gsl_vector * x, void *params) {
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf f;
    int             status;
    int             iter = 0;
    int             p = x->size;

    f.f = Func_f;
    f.df = Func_df;
    f.fdf = Func_fdf;
    f.n = p;
    f.params = params;

    T = gsl_multimin_fdfminimizer_steepest_descent;
    s = gsl_multimin_fdfminimizer_alloc(T, p);

    gsl_multimin_fdfminimizer_set(s, &f, x, 0.1, 0.1);

    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status)
            break;

        if (printFlag && iter % 10 == 0)
            printf("SteepDescMin: %3d    %8.3lf     %8.3lf\n", iter, s->f, gsl_blas_dnrm2(s->gradient));

        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
    } while (status == GSL_CONTINUE && iter < 100);

    if (printFlag) {
        printf("SteepDescMin: %3d    %8.3lf     %8.3lf\n", iter, s->f, gsl_blas_dnrm2(s->gradient));
        printf("\n");
    }

    gsl_vector_memcpy(x, s->x);

    gsl_multimin_fdfminimizer_free(s);

    return s->f;
}

double SimplexMin(double (*Func_f) (const gsl_vector *, void *), gsl_vector * x, void *params) {
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function f;
    gsl_vector     *step;
    int             status;
    double          size = 0.0;
    int             iter = 0;
    int             p = x->size;

    step = gsl_vector_alloc(p);
    gsl_vector_set_all(step, 1e-1);

    f.f = Func_f;
    f.n = p;
    f.params = params;

    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc(T, p);

    gsl_multimin_fminimizer_set(s, &f, x, step);

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-5);

        if (printFlag && iter % 1000 == 0)
            printf("LMIN: %3d    %8.3lf     %8.3lf\n", iter, s->fval, size);

    } while (status == GSL_CONTINUE && iter < 50000);

    if (printFlag)
        printf("LMIN: %3d    %8.3lf     %8.3lf\n", iter, s->fval, size);

    gsl_vector_memcpy(x, s->x);

    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(s);

    return s->fval;
}

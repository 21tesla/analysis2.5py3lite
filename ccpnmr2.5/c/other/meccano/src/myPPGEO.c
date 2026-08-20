#include "../inc/myPPGEO.h"

/* Peptide Plane */
static double   _ca1[] = { 0.000, 0.000, 0.000 };
static double   _c1_[] = { 1.430, -0.529, 0.000 };
static double   _o1_[] = { 1.669, -1.746, 0.000 };
static double   _n2_[] = { 2.404, 0.382, 0.000 };
static double   _hn2[] = { 2.161, 1.373, 0.000 };
static double   _ca2[] = { 3.806, -0.006, 0.000 };

pPlane_t       *InitPepPla() {
    gsl_vector     *tmp = gsl_vector_calloc(3);

    pPlane_t       *pPlane = (pPlane_t *) malloc(sizeof(pPlane_t));

    gsl_vector_view v_ca1_ref = gsl_vector_view_array(_ca1, 3);
    gsl_vector_view v_c1__ref = gsl_vector_view_array(_c1_, 3);
    gsl_vector_view v_n2__ref = gsl_vector_view_array(_n2_, 3);
    gsl_vector_view v_hn2_ref = gsl_vector_view_array(_hn2, 3);
    gsl_vector_view v_ca2_ref = gsl_vector_view_array(_ca2, 3);
    gsl_vector_view v_o1__ref = gsl_vector_view_array(_o1_, 3);

    pPlane->ca1_ref = gsl_vector_calloc(3);
    pPlane->c1__ref = gsl_vector_calloc(3);
    pPlane->n2__ref = gsl_vector_calloc(3);
    pPlane->hn2_ref = gsl_vector_calloc(3);
    pPlane->ca2_ref = gsl_vector_calloc(3);
    pPlane->o1__ref = gsl_vector_calloc(3);

    gsl_vector_memcpy(pPlane->ca1_ref, &v_ca1_ref.vector);
    gsl_vector_memcpy(pPlane->c1__ref, &v_c1__ref.vector);
    gsl_vector_memcpy(pPlane->n2__ref, &v_n2__ref.vector);
    gsl_vector_memcpy(pPlane->hn2_ref, &v_hn2_ref.vector);
    gsl_vector_memcpy(pPlane->ca2_ref, &v_ca2_ref.vector);
    gsl_vector_memcpy(pPlane->o1__ref, &v_o1__ref.vector);

    pPlane->c_n__ref = gsl_vector_calloc(3);
    pPlane->n_ca_ref = gsl_vector_calloc(3);
    pPlane->ca_c_ref = gsl_vector_calloc(3);
    Points2UVector(pPlane->c1__ref, pPlane->n2__ref, pPlane->c_n__ref);
    Points2UVector(pPlane->n2__ref, pPlane->ca2_ref, pPlane->n_ca_ref);
    Points2UVector(pPlane->ca1_ref, pPlane->c1__ref, pPlane->ca_c_ref);

    pPlane->x_pp_ref = gsl_vector_calloc(3);
    pPlane->y_pp_ref = gsl_vector_calloc(3);
    pPlane->z_pp_ref = gsl_vector_calloc(3);

    Points2UVector(pPlane->ca1_ref, pPlane->ca2_ref, pPlane->z_pp_ref);
    Points2UVector(pPlane->c1__ref, pPlane->o1__ref, tmp);
    CrossProduct(pPlane->z_pp_ref, tmp, pPlane->y_pp_ref);
    Normalize(pPlane->y_pp_ref);
    CrossProduct(pPlane->y_pp_ref, pPlane->z_pp_ref, pPlane->x_pp_ref);
    Normalize(pPlane->x_pp_ref);

    gsl_vector_free(tmp);

    printf("-------------------------------------------------------\n");
    printf("%10s %10s %10s\n", "Distances", "In PPGEO", "In cstRDC");
    printf("%10s %10.3lf %10s\n", "CB", R_CB, "X");
    printf("%10s %10.3lf %10s\n", "HA", R_HA, "X");
    printf("%10s %10.3lf %10.3lf\n", "CN", R_CN, Dist(pPlane->c1__ref, pPlane->n2__ref));
    printf("%10s %10.3lf %10.3lf\n", "CC", R_CC, Dist(pPlane->c1__ref, pPlane->ca1_ref));
    printf("%10s %10.3lf %10.3lf\n", "CH", R_CH, Dist(pPlane->c1__ref, pPlane->hn2_ref));
    printf("%10s %10.3lf %10.3lf\n", "NH", R_NH, Dist(pPlane->n2__ref, pPlane->hn2_ref));
    printf("%10s %10.3lf %10.3lf\n", "CA1H", R_CA1H, Dist(pPlane->ca1_ref, pPlane->hn2_ref));
    printf("%10s %10.3lf %10.3lf\n", "CA2H", R_CA2H, Dist(pPlane->ca2_ref, pPlane->hn2_ref));
    printf("-------------------------------------------------------\n");

    return pPlane;
}

pPlane_t       *InitPepPlaFwd() {
    pPlane_t       *pPlane = InitPepPla();
    gsl_vector_view origin = gsl_vector_view_array(_ca1, 3);

    gsl_vector_sub(pPlane->ca1_ref, &origin.vector);
    gsl_vector_sub(pPlane->c1__ref, &origin.vector);
    gsl_vector_sub(pPlane->n2__ref, &origin.vector);
    gsl_vector_sub(pPlane->hn2_ref, &origin.vector);
    gsl_vector_sub(pPlane->ca2_ref, &origin.vector);
    gsl_vector_sub(pPlane->o1__ref, &origin.vector);

    return pPlane;
}

pPlane_t       *InitPepPlaRev() {
    pPlane_t       *pPlane = InitPepPla();
    gsl_vector_view origin = gsl_vector_view_array(_ca2, 3);

    gsl_vector_sub(pPlane->ca1_ref, &origin.vector);
    gsl_vector_sub(pPlane->c1__ref, &origin.vector);
    gsl_vector_sub(pPlane->n2__ref, &origin.vector);
    gsl_vector_sub(pPlane->hn2_ref, &origin.vector);
    gsl_vector_sub(pPlane->ca2_ref, &origin.vector);
    gsl_vector_sub(pPlane->o1__ref, &origin.vector);

    return pPlane;
}

void DelPepPla(pPlane_t *pPlane)
{
    gsl_vector_free(pPlane->ca1_ref);
    gsl_vector_free(pPlane->c1__ref);
    gsl_vector_free(pPlane->n2__ref);
    gsl_vector_free(pPlane->hn2_ref);
    gsl_vector_free(pPlane->ca2_ref);
    gsl_vector_free(pPlane->o1__ref);
    gsl_vector_free(pPlane->c_n__ref);
    gsl_vector_free(pPlane->n_ca_ref);
    gsl_vector_free(pPlane->ca_c_ref);
    gsl_vector_free(pPlane->x_pp_ref);
    gsl_vector_free(pPlane->y_pp_ref);
    gsl_vector_free(pPlane->z_pp_ref);

    free(pPlane);
}


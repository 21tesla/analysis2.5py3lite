#include "../inc/myRAMACHANDRAN.h"

extern gsl_rng *r;

// ramachandran_t *InitRama(int resolution) {
//     ramachandran_t *rama = (ramachandran_t *) malloc(sizeof(ramachandran_t));
// 
//     rama->resolution = resolution;
//     rama->gen = gsl_matrix_alloc(resolution, resolution);
//     rama->prp = gsl_matrix_alloc(resolution, resolution);
//     rama->pro = gsl_matrix_alloc(resolution, resolution);
//     rama->gly = gsl_matrix_alloc(resolution, resolution);
// 
//     return rama;
// }

// void ReadRamaFile(char *fileName, gsl_matrix * ramaMat, int resolution) {
//     int             i, j;
//     double          ph, ps, h;
//     FILE           *fileIn;
// 
//     if (!(fileIn = fopen(fileName, "r"))) {
//         printf("Can't open %s", fileName);
//         exit(0);
//     }
//     while (!feof(fileIn)) {
//         fscanf(fileIn, "%lf %lf %lf\n", &ph, &ps, &h);
//         i = (int) ((ph + 180.) / 360. * resolution);
//         j = (int) ((ps + 180.) / 360. * resolution);
//         gsl_matrix_set(ramaMat, i, j, h);
//     }
//     fclose(fileIn);
// }

// void ReadRamaFiles(ramachandran_t * rama) {
//     int             i, j;
//     double          n_gen = 0., n_gly = 0., n_pro = 0., n_prp = 0.;
// 
//     ReadRamaFile("gen.hist", rama->gen, rama->resolution);
//     ReadRamaFile("prp.hist", rama->prp, rama->resolution);
//     ReadRamaFile("pro.hist", rama->pro, rama->resolution);
//     ReadRamaFile("gly.hist", rama->gly, rama->resolution);
// 
//     for (i = 0; i < rama->resolution; i++)
//         for (j = 0; j < rama->resolution; j++) {
//             n_gen += gsl_matrix_get(rama->gen, i, j);
//             n_gly += gsl_matrix_get(rama->gly, i, j);
//             n_pro += gsl_matrix_get(rama->pro, i, j);
//             n_prp += gsl_matrix_get(rama->prp, i, j);
//         }
//     gsl_matrix_scale(rama->gen, 1. / n_gen);
//     gsl_matrix_scale(rama->gly, 1. / n_gly);
//     gsl_matrix_scale(rama->pro, 1. / n_pro);
//     gsl_matrix_scale(rama->prp, 1. / n_prp);
// 
//     for (i = 0; i < rama->resolution; i++)
//         for (j = 0; j < rama->resolution; j++) {
//             gsl_matrix_set(rama->gen, i, j, gsl_matrix_get(rama->gen, i, j) ? -log(gsl_matrix_get(rama->gen, i, j)) : 1e16);
//             gsl_matrix_set(rama->gly, i, j, gsl_matrix_get(rama->gly, i, j) ? -log(gsl_matrix_get(rama->gly, i, j)) : 1e16);
//             gsl_matrix_set(rama->pro, i, j, gsl_matrix_get(rama->pro, i, j) ? -log(gsl_matrix_get(rama->pro, i, j)) : 1e16);
//             gsl_matrix_set(rama->prp, i, j, gsl_matrix_get(rama->prp, i, j) ? -log(gsl_matrix_get(rama->prp, i, j)) : 1e16);
//         }
// 
// }

// double ERamachandran(ramachandran_t * rama, polyPeptide_t * polyP, int ii) {
//     int             i, j;
//     double          energy;
//     double          phiSh = polyP->acAm[ii].phi + PI, psiSh = polyP->acAm[ii].psi + PI;
// 
//     gsl_sf_angle_restrict_pos_e(&phiSh);
//     gsl_sf_angle_restrict_pos_e(&psiSh);
// 
//     i = (int) (phiSh / (2. * PI) * rama->resolution);
//     j = (int) (psiSh / (2. * PI) * rama->resolution);
// 
//     if (i == 180)
//         i = 0;
//     if (j == 180)
//         j = 0;
// 
//     if (polyP->acAm[ii + 1].acAmType == PRO)
//         energy = gsl_matrix_get(rama->prp, i, j);
//     else if (polyP->acAm[ii].acAmType == GLY)
//         energy = gsl_matrix_get(rama->gly, i, j);
//     else if (polyP->acAm[ii].acAmType == PRO)
//         energy = gsl_matrix_get(rama->pro, i, j);
//     else
//         energy = gsl_matrix_get(rama->gen, i, j);
// 
//     return energy;
// }

ramaDB_t       *ReadRamaDBFile(char *ramaDBFileName) {
    char            resName[4];
    double          phi_p, psi_p;
    acAmType_t      acAmType;
    FILE           *fileIn;
    ramaDB_t       *ramaDB = (ramaDB_t *) malloc(sizeof(ramaDB_t));

    for (acAmType = 0; acAmType < NB_AC_AM_TYPE; acAmType++)
        ramaDB->pointNb[acAmType] = 0;

    if (!(fileIn = fopen(ramaDBFileName, "r"))) {
        printf("Can't open %s", ramaDBFileName);
        return NULL;
    }
    while (!feof(fileIn)) {
        fscanf(fileIn, "%s %lf %lf\n", resName, &phi_p, &psi_p);
        AcAmName2Type(resName, &acAmType);
        if (acAmType < NB_AC_AM_TYPE)
            ramaDB->pointNb[acAmType]++;
    }
    fclose(fileIn);

    for (acAmType = 0; acAmType < NB_AC_AM_TYPE; acAmType++) {
        ramaDB->phi[acAmType] = gsl_vector_alloc(ramaDB->pointNb[acAmType]);
        ramaDB->psi[acAmType] = gsl_vector_alloc(ramaDB->pointNb[acAmType]);
    }

    for (acAmType = 0; acAmType < NB_AC_AM_TYPE; acAmType++)
        ramaDB->pointNb[acAmType] = 0;

    fileIn = fopen(ramaDBFileName, "r");
    while (!feof(fileIn)) {
        fscanf(fileIn, "%s %lf %lf\n", resName, &phi_p, &psi_p);
        AcAmName2Type(resName, &acAmType);
        if (acAmType < NB_AC_AM_TYPE) {
            gsl_vector_set(ramaDB->phi[acAmType], ramaDB->pointNb[acAmType], DEG2RAD(phi_p));
            gsl_vector_set(ramaDB->psi[acAmType], ramaDB->pointNb[acAmType], DEG2RAD(psi_p));
            ramaDB->pointNb[acAmType]++;
        }
    }
    fclose(fileIn);

    return ramaDB;
}

void AleaRamaDB(ramaDB_t * RamaDB, polyPeptide_t * polyP, int ii, double *phi, double *psi) {
    int             i;
    acAmType_t      acAmType;

    if (polyP->acAm[ii + 1].acAmType == PRO) {
        if (polyP->acAm[ii].acAmType == PRO)
            acAmType = PPR;
        else
            acAmType = PRP;
    } else
        acAmType = polyP->acAm[ii].acAmType;

    i = gsl_rng_uniform_int(r, RamaDB->pointNb[acAmType]);
    *phi = gsl_vector_get(RamaDB->phi[acAmType], i);
    *psi = gsl_vector_get(RamaDB->psi[acAmType], i);
}

void DelRamaDB(ramaDB_t * RamaDB)
{
    acAmType_t      acAmType;

    for (acAmType = 0; acAmType < NB_AC_AM_TYPE; acAmType++) {
        gsl_vector_free(RamaDB->phi[acAmType]);
        gsl_vector_free(RamaDB->psi[acAmType]);
    }

    free(RamaDB);
}


#include "../inc/myJUNC.h"

double CalcChi2_JuncGeo(polyPeptide_t * polyP, data_t * data, int resI) {
    double          chi2 = 0.0;

    if (data->pepPla[resI].phiFlag)
        chi2 += CHI2/*_FLAT*/(RAD2DEG(polyP->acAm[resI].phi), data->pepPla[resI].phi0, data->pepPla[resI].phiEr);

    if (data->pepPla[resI].psiFlag)
        chi2 += CHI2/*_FLAT*/(RAD2DEG(polyP->acAm[resI].psi), data->pepPla[resI].psi0, data->pepPla[resI].psiEr);

    if (polyP->tetFlag)
        chi2 += CHI2(polyP->acAm[resI].tet, polyP->tet0, polyP->tetEr);

    return chi2;
}

void DisplayJunctionParam(FILE * output, polyPeptide_t * polyP, int resI) {
    char            resName[4];

    AcAmType2Name(polyP->acAm[resI].acAmType, resName);
    fprintf(output, "\nJunction Parameters:\n");
    fprintf(output, "    %s %d  /  Phi: %6.1lf   Psi: %6.1lf   Tet: %6.1lf\n\n", resName, polyP->acAm[resI].resNb, RAD2DEG(polyP->acAm[resI].phi), RAD2DEG(polyP->acAm[resI].psi), RAD2DEG(polyP->acAm[resI].tet));
}

double DisplayChi2_JuncGeo(FILE * output, polyPeptide_t * polyP, data_t * data, int resI) {
    double          chi2 = 0.0;
    double          chi2Phi = 0.0, chi2Psi = 0.0, chi2Tet = 0.0;

    fprintf(output, "%-15s %6s %6s %6s\n", "Junc Geo Chi2:", "Phi", "Psi", "Tet");

    if (data->pepPla[resI].phiFlag) {
        chi2 += chi2Psi = CHI2(RAD2DEG(polyP->acAm[resI].phi), data->pepPla[resI].phi0, data->pepPla[resI].phiEr);
        fprintf(output, "%15s %6.1lf", "", chi2Psi);
    } else
        fprintf(output, "%15s %6s", "", "XXX");

    if (data->pepPla[resI].psiFlag) {
        chi2 += chi2Phi = CHI2(RAD2DEG(polyP->acAm[resI].psi), data->pepPla[resI].psi0, data->pepPla[resI].psiEr);
        fprintf(output, " %6.1lf", chi2Phi);
    } else
        fprintf(output, " %6s", "XXX");

    if (polyP->tetFlag) {
        chi2 += chi2Tet = CHI2(polyP->acAm[resI].tet, polyP->tet0, polyP->tetEr);
        fprintf(output, " %6.1lf", chi2Tet);
    } else
        fprintf(output, " %6s", "XXX");

    fprintf(output, "\n\n");
    return chi2;
}

#include "../inc/mySTRUCT.h"

extern int      printFlag;

polyPeptide_t  *InitPolyPeptide(int totalPepPlaNb, int firstPepPlaNb) {
    int             i;
    polyPeptide_t  *polyP = (polyPeptide_t *) malloc(sizeof(polyPeptide_t));

    polyP->totalPepPlaNb = totalPepPlaNb;
    polyP->firstPepPlaNb = firstPepPlaNb;

    polyP->tet0 = 111.0;
    polyP->dTet0 = 4.0;
    polyP->tetEr = 4.0;
    polyP->tetFlag = 1;

    polyP->acAm = (acAm_t *) malloc((polyP->totalPepPlaNb + 2) * sizeof(acAm_t));

    for (i = 0; i <= polyP->totalPepPlaNb; i++) {
        polyP->acAm[i].resNb = i + firstPepPlaNb - 1;

        polyP->acAm[i].n_ = gsl_vector_alloc(3);
        polyP->acAm[i].h_ = gsl_vector_alloc(3);
        polyP->acAm[i].ca = gsl_vector_alloc(3);
        polyP->acAm[i].ha = gsl_vector_alloc(3);
        polyP->acAm[i].cb = gsl_vector_alloc(3);
        polyP->acAm[i].c_ = gsl_vector_alloc(3);
        polyP->acAm[i].o_ = gsl_vector_alloc(3);

        polyP->acAm[i].xAxis = gsl_vector_alloc(3);
        polyP->acAm[i].yAxis = gsl_vector_alloc(3);
        polyP->acAm[i].zAxis = gsl_vector_alloc(3);

        polyP->acAm[i].n_Flag = 0;
        polyP->acAm[i].h_Flag = 0;
        polyP->acAm[i].caFlag = 0;
        polyP->acAm[i].haFlag = 0;
        polyP->acAm[i].cbFlag = 0;
        polyP->acAm[i].c_Flag = 0;
        polyP->acAm[i].o_Flag = 0;
    }
    polyP->acAm[polyP->totalPepPlaNb + 1].resNb = -1;

    return polyP;
}

void DelPolyPeptide(polyPeptide_t * polyP) {
    int             i;

    for (i = 0; i <= polyP->totalPepPlaNb; i++) {
        gsl_vector_free(polyP->acAm[i].n_);
        gsl_vector_free(polyP->acAm[i].ca);
        gsl_vector_free(polyP->acAm[i].c_);
        gsl_vector_free(polyP->acAm[i].h_);
        gsl_vector_free(polyP->acAm[i].o_);
        gsl_vector_free(polyP->acAm[i].cb);
        gsl_vector_free(polyP->acAm[i].ha);

        gsl_vector_free(polyP->acAm[i].xAxis);
        gsl_vector_free(polyP->acAm[i].yAxis);
        gsl_vector_free(polyP->acAm[i].zAxis);
    }
    free(polyP->acAm);
    free(polyP);
}

void AddCbHa(int ppn, polyPeptide_t * polyP) {
    double          q1 = DEG2RAD(-125.25);

    gsl_vector     *ca_n = gsl_vector_alloc(3);
    gsl_vector     *ca_c = gsl_vector_alloc(3);
    gsl_vector     *perp = gsl_vector_alloc(3);
    gsl_vector     *para = gsl_vector_alloc(3);
    gsl_vector     *cbca = gsl_vector_alloc(3);
    gsl_vector     *haca = gsl_vector_alloc(3);

    Points2Vector(polyP->acAm[ppn].ca, polyP->acAm[ppn].n_, ca_n);
    Points2Vector(polyP->acAm[ppn].ca, polyP->acAm[ppn].c_, ca_c);

    /* 
     * caLCULATE RDC's OF ca-cb   
     */

    CrossProduct(ca_c, ca_n, perp);
    Normalize(perp);

    gsl_vector_memcpy(para, ca_n);
    gsl_vector_add(para, ca_c);
    Normalize(para);

    gsl_vector_scale(para, cos(q1));
    gsl_vector_scale(perp, sin(q1));

    gsl_vector_memcpy(cbca, para);
    gsl_vector_add(cbca, perp);
    Normalize(cbca);
    ScalingTranslation(R_CB, polyP->acAm[ppn].ca, cbca, polyP->acAm[ppn].cb);

    gsl_vector_memcpy(haca, para);
    gsl_vector_sub(haca, perp);
    Normalize(haca);
    ScalingTranslation(R_HA, polyP->acAm[ppn].ca, haca, polyP->acAm[ppn].ha);

    polyP->acAm[ppn].cbFlag = polyP->acAm[ppn].haFlag = 1;

    gsl_vector_free(ca_n);
    gsl_vector_free(ca_c);
    gsl_vector_free(perp);
    gsl_vector_free(para);
    gsl_vector_free(cbca);
    gsl_vector_free(haca);
}

void OrientatePepPlaFwd(const gsl_vector * x, int firstPepPlaFlag, int ppn, pPlane_t * pPlane, polyPeptide_t * polyP) {
    gsl_matrix     *rotMat = gsl_matrix_calloc(3, 3);

    if (firstPepPlaFlag) {
        AxisAngle2Mat(x, rotMat);
    } else {
        double          phi, psi, tet;
        double          phi0, psi0, tet0;

        gsl_matrix     *rotPhi = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rotPsi = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rotTet = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rottmp = gsl_matrix_alloc(3, 3);

        gsl_vector     *c_n__pre = gsl_vector_calloc(3);
        gsl_vector     *n_ca_pre = gsl_vector_calloc(3);
        gsl_vector     *ca_c_tmp = gsl_vector_calloc(3);
        gsl_vector     *ca_c = gsl_vector_calloc(3);
        gsl_vector     *normal = gsl_vector_calloc(3);

        phi = GET(x, 0);
        psi = GET(x, 1);
        tet = GET(x, 2);

        Points2Vector(polyP->acAm[ppn - 2].c_, polyP->acAm[ppn - 1].n_, c_n__pre);
        Points2Vector(polyP->acAm[ppn - 1].n_, polyP->acAm[ppn - 1].ca, n_ca_pre);

        CalcPhiPsiTet(c_n__pre, n_ca_pre, pPlane->ca_c_ref, pPlane->c_n__ref, &phi0, &psi0, &tet0);

        CrossProduct(pPlane->ca_c_ref, n_ca_pre, normal);
        RotMat(normal, -(tet0 - tet), rotTet);
        RotMat(n_ca_pre, -(phi0 - phi), rotPhi);
        Rotation(rotTet, pPlane->ca_c_ref, ca_c_tmp);
        Rotation(rotPhi, ca_c_tmp, ca_c);
        RotMat(ca_c, -(psi0 - psi), rotPsi);

        MatrixMultiplication(rotPhi, rotTet, rottmp);
        MatrixMultiplication(rotPsi, rottmp, rotMat);

        polyP->acAm[ppn - 1].phi = phi;
        polyP->acAm[ppn - 1].psi = psi;
        polyP->acAm[ppn - 1].tet = tet;

        gsl_matrix_free(rotPhi);
        gsl_matrix_free(rotPsi);
        gsl_matrix_free(rotTet);
        gsl_matrix_free(rottmp);

        gsl_vector_free(c_n__pre);
        gsl_vector_free(n_ca_pre);
        gsl_vector_free(ca_c_tmp);
        gsl_vector_free(ca_c);
        gsl_vector_free(normal);
    }

    RotationTranslation(rotMat, polyP->acAm[ppn - 1].ca, pPlane->c1__ref, polyP->acAm[ppn - 1].c_);
    RotationTranslation(rotMat, polyP->acAm[ppn - 1].ca, pPlane->n2__ref, polyP->acAm[ppn].n_);
    RotationTranslation(rotMat, polyP->acAm[ppn - 1].ca, pPlane->hn2_ref, polyP->acAm[ppn].h_);
    RotationTranslation(rotMat, polyP->acAm[ppn - 1].ca, pPlane->o1__ref, polyP->acAm[ppn - 1].o_);
    RotationTranslation(rotMat, polyP->acAm[ppn - 1].ca, pPlane->ca2_ref, polyP->acAm[ppn].ca);

    polyP->acAm[ppn - 1].caFlag = polyP->acAm[ppn - 1].c_Flag = polyP->acAm[ppn].n_Flag = polyP->acAm[ppn].h_Flag = polyP->acAm[ppn - 1].o_Flag = polyP->acAm[ppn].caFlag = 1;

    if (!firstPepPlaFlag) {
        AddCbHa(ppn - 1, polyP);
    }

    Rotation(rotMat, pPlane->x_pp_ref, polyP->acAm[ppn].xAxis);
    Rotation(rotMat, pPlane->y_pp_ref, polyP->acAm[ppn].yAxis);
    Rotation(rotMat, pPlane->z_pp_ref, polyP->acAm[ppn].zAxis);

    gsl_matrix_free(rotMat);
}

void OrientatePepPlaRev(const gsl_vector * x, int firstPepPlaFlag, int ppn, pPlane_t * pPlane, polyPeptide_t * polyP) {
    gsl_matrix     *rotMat = gsl_matrix_calloc(3, 3);

    if (firstPepPlaFlag) {
        AxisAngle2Mat(x, rotMat);
    } else {
        double          phi, psi, tet;
        double          phi0, psi0, tet0;

        gsl_matrix     *rotPhi = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rotPsi = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rotTet = gsl_matrix_alloc(3, 3);
        gsl_matrix     *rottmp = gsl_matrix_alloc(3, 3);

        gsl_vector     *ca_c_pre = gsl_vector_calloc(3);
        gsl_vector     *c_n__pre = gsl_vector_calloc(3);
        gsl_vector     *n_ca_tmp = gsl_vector_calloc(3);
        gsl_vector     *n_ca = gsl_vector_calloc(3);
        gsl_vector     *normal = gsl_vector_calloc(3);

        phi = GET(x, 0);
        psi = GET(x, 1);
        tet = GET(x, 2);

        Points2Vector(polyP->acAm[ppn].ca, polyP->acAm[ppn].c_, ca_c_pre);
        Points2Vector(polyP->acAm[ppn].c_, polyP->acAm[ppn + 1].n_, c_n__pre);

        CalcPhiPsiTet(pPlane->c_n__ref, pPlane->n_ca_ref, ca_c_pre, c_n__pre, &phi0, &psi0, &tet0);

        CrossProduct(ca_c_pre, pPlane->n_ca_ref, normal);
        RotMat(normal, (tet0 - tet), rotTet);
        RotMat(ca_c_pre, (psi0 - psi), rotPsi);
        Rotation(rotTet, pPlane->n_ca_ref, n_ca_tmp);
        Rotation(rotPsi, n_ca_tmp, n_ca);
        RotMat(n_ca, (phi0 - phi), rotPhi);

        MatrixMultiplication(rotPsi, rotTet, rottmp);
        MatrixMultiplication(rotPhi, rottmp, rotMat);

        polyP->acAm[ppn].phi = phi;
        polyP->acAm[ppn].psi = psi;
        polyP->acAm[ppn].tet = tet;

        gsl_matrix_free(rotPhi);
        gsl_matrix_free(rotPsi);
        gsl_matrix_free(rotTet);
        gsl_matrix_free(rottmp);

        gsl_vector_free(ca_c_pre);
        gsl_vector_free(c_n__pre);
        gsl_vector_free(n_ca_tmp);
        gsl_vector_free(n_ca);
        gsl_vector_free(normal);
    }

    RotationTranslation(rotMat, polyP->acAm[ppn].ca, pPlane->ca1_ref, polyP->acAm[ppn - 1].ca);
    RotationTranslation(rotMat, polyP->acAm[ppn].ca, pPlane->c1__ref, polyP->acAm[ppn - 1].c_);
    RotationTranslation(rotMat, polyP->acAm[ppn].ca, pPlane->n2__ref, polyP->acAm[ppn].n_);
    RotationTranslation(rotMat, polyP->acAm[ppn].ca, pPlane->hn2_ref, polyP->acAm[ppn].h_);
    RotationTranslation(rotMat, polyP->acAm[ppn].ca, pPlane->o1__ref, polyP->acAm[ppn - 1].o_);

    polyP->acAm[ppn].caFlag = polyP->acAm[ppn - 1].caFlag = polyP->acAm[ppn - 1].c_Flag = polyP->acAm[ppn].n_Flag = polyP->acAm[ppn].h_Flag = polyP->acAm[ppn - 1].o_Flag = 1;

    if (!firstPepPlaFlag) {
        AddCbHa(ppn, polyP);
    }

    Rotation(rotMat, pPlane->x_pp_ref, polyP->acAm[ppn].xAxis);
    Rotation(rotMat, pPlane->y_pp_ref, polyP->acAm[ppn].yAxis);
    Rotation(rotMat, pPlane->z_pp_ref, polyP->acAm[ppn].zAxis);

    gsl_matrix_free(rotMat);
}

void ReadSequenceFile(char *fileName, polyPeptide_t * polyP) {
    int             seq_num;
    char            seq_name[2];
    FILE           *fileIn;

    if (!(fileIn = fopen(fileName, "r"))) {
        printf("Can't open: %s", fileName);
        exit(0);
    }

    while (!feof(fileIn)) {

        fscanf(fileIn, "%d %s\n", &seq_num, seq_name);
        SetSequence(polyP, seq_num, seq_name);
    }

    fclose(fileIn);
}

void SetSequence(polyPeptide_t * polyP, int seq_num, char *seq_name) {
    int             i;

    i = 0;
    while (polyP->acAm[i].resNb != seq_num && polyP->acAm[i].resNb != -1)
        i++;

    if (!strcmp(seq_name, "G"))
        polyP->acAm[i].acAmType = GLY;
    else if (!strcmp(seq_name, "A"))
        polyP->acAm[i].acAmType = ALA;
    else if (!strcmp(seq_name, "R"))
        polyP->acAm[i].acAmType = ARG;
    else if (!strcmp(seq_name, "N"))
        polyP->acAm[i].acAmType = ASN;
    else if (!strcmp(seq_name, "D"))
        polyP->acAm[i].acAmType = ASP;
    else if (!strcmp(seq_name, "C"))
        polyP->acAm[i].acAmType = CYS;
    else if (!strcmp(seq_name, "Q"))
        polyP->acAm[i].acAmType = GLN;
    else if (!strcmp(seq_name, "E"))
        polyP->acAm[i].acAmType = GLU;
    else if (!strcmp(seq_name, "H"))
        polyP->acAm[i].acAmType = HIS;
    else if (!strcmp(seq_name, "I"))
        polyP->acAm[i].acAmType = ILE;
    else if (!strcmp(seq_name, "L"))
        polyP->acAm[i].acAmType = LEU;
    else if (!strcmp(seq_name, "K"))
        polyP->acAm[i].acAmType = LYS;
    else if (!strcmp(seq_name, "M"))
        polyP->acAm[i].acAmType = MET;
    else if (!strcmp(seq_name, "F"))
        polyP->acAm[i].acAmType = PHE;
    else if (!strcmp(seq_name, "P"))
        polyP->acAm[i].acAmType = PRO;
    else if (!strcmp(seq_name, "S"))
        polyP->acAm[i].acAmType = SER;
    else if (!strcmp(seq_name, "T"))
        polyP->acAm[i].acAmType = THR;
    else if (!strcmp(seq_name, "W"))
        polyP->acAm[i].acAmType = TRP;
    else if (!strcmp(seq_name, "Y"))
        polyP->acAm[i].acAmType = TYR;
    else if (!strcmp(seq_name, "V"))
        polyP->acAm[i].acAmType = VAL;
}

void AcAmType2Name(acAmType_t acAmType, char *name) {
    switch (acAmType) {
        case GLY:
            strcpy(name, "GLY");
            break;
        case ALA:
            strcpy(name, "ALA");
            break;
        case ARG:
            strcpy(name, "ARG");
            break;
        case ASN:
            strcpy(name, "ASN");
            break;
        case ASP:
            strcpy(name, "ASP");
            break;
        case CYS:
            strcpy(name, "CYS");
            break;
        case GLN:
            strcpy(name, "GLN");
            break;
        case GLU:
            strcpy(name, "GLU");
            break;
        case HIS:
            strcpy(name, "HIS");
            break;
        case ILE:
            strcpy(name, "ILE");
            break;
        case LEU:
            strcpy(name, "LEU");
            break;
        case LYS:
            strcpy(name, "LYS");
            break;
        case MET:
            strcpy(name, "MET");
            break;
        case PHE:
            strcpy(name, "PHE");
            break;
        case SER:
            strcpy(name, "SER");
            break;
        case THR:
            strcpy(name, "THR");
            break;
        case TRP:
            strcpy(name, "TRP");
            break;
        case TYR:
            strcpy(name, "TYR");
            break;
        case VAL:
            strcpy(name, "VAL");
            break;
        case PRO:
            strcpy(name, "PRO");
            break;
        case PRP:
            strcpy(name, "PRP");
            break;
        case PPR:
            strcpy(name, "PPR");
            break;
       default:
            break;
    }
}

void AcAmName2Type(char name[4], acAmType_t * acAmType) {
    if (!strcmp(name, "GLY"))
        (*acAmType) = GLY;
    else if (!strcmp(name, "ALA"))
        (*acAmType) = ALA;
    else if (!strcmp(name, "ARG"))
        (*acAmType) = ARG;
    else if (!strcmp(name, "ASN"))
        (*acAmType) = ASN;
    else if (!strcmp(name, "ASP"))
        (*acAmType) = ASP;
    else if (!strcmp(name, "CYS"))
        (*acAmType) = CYS;
    else if (!strcmp(name, "GLN"))
        (*acAmType) = GLN;
    else if (!strcmp(name, "GLU"))
        (*acAmType) = GLU;
    else if (!strcmp(name, "HIS"))
        (*acAmType) = HIS;
    else if (!strcmp(name, "ILE"))
        (*acAmType) = ILE;
    else if (!strcmp(name, "LEU"))
        (*acAmType) = LEU;
    else if (!strcmp(name, "LYS"))
        (*acAmType) = LYS;
    else if (!strcmp(name, "MET"))
        (*acAmType) = MET;
    else if (!strcmp(name, "PHE"))
        (*acAmType) = PHE;
    else if (!strcmp(name, "SER"))
        (*acAmType) = SER;
    else if (!strcmp(name, "THR"))
        (*acAmType) = THR;
    else if (!strcmp(name, "TRP"))
        (*acAmType) = TRP;
    else if (!strcmp(name, "TYR"))
        (*acAmType) = TYR;
    else if (!strcmp(name, "VAL"))
        (*acAmType) = VAL;
    else if (!strcmp(name, "PRO"))
        (*acAmType) = PRO;
    else if (!strcmp(name, "PRP"))
        (*acAmType) = PRP;
    else if (!strcmp(name, "PPR"))
        (*acAmType) = PPR;
    else
        (*acAmType) = NB_AC_AM_TYPE + 1;
}

polyPeptide_t  *ReadPdbFile(char *fileName) {
    polyPeptide_t  *polyP;
    FILE           *fileIn;

    int             atn, resn;
    double          x, y, z;
    char            at[7], att[5], res_name[4], end[100];

    int             firstPepPlaNb = 1000000;
    int             Last_pPlane = 0;
    int             totalPepPlaNb = 0;

//     printf("Read: %s\n", fileName);
    if (!(fileIn = fopen(fileName, "r"))) {
        printf("Can't open: %s", fileName);
        exit(0);
    }

    while (!feof(fileIn)) {
        fscanf(fileIn, "%s %d %s %s %d %lf %lf %lf %[^\n] \n", at, &atn, att, res_name, &resn, &x, &y, &z, end);

        if ((!strcmp(at, "ATOM")) && (!strcmp(att, "N"))) {
            firstPepPlaNb = GSL_MIN(firstPepPlaNb, resn);
            Last_pPlane = GSL_MAX(Last_pPlane, resn);
        }
    }
    fclose(fileIn);
    fileIn = NULL;

    totalPepPlaNb = Last_pPlane - firstPepPlaNb + 1;

//     fprintf(stderr, "firstPepPlaNb: %d   totalPepPlaNb: %d\n", firstPepPlaNb, totalPepPlaNb);

    polyP = InitPolyPeptide(totalPepPlaNb, firstPepPlaNb);
    polyP->windowOffset = 0;
    polyP->windowSize = totalPepPlaNb;

    if (!(fileIn = fopen(fileName, "r"))) {
        printf("Can't open: %s", fileName);
        exit(0);
    }
    while (!feof(fileIn)) {
        acAmType_t      Res_Type;
        int             offset = firstPepPlaNb - 1;

        fscanf(fileIn, "%s %d %s %s %d %lf %lf %lf %[^\n] \n", at, &atn, att, res_name, &resn, &x, &y, &z, end);

        if (!strcmp(att, "CA")) {
            polyP->acAm[resn - offset].resNb = resn;
            AcAmName2Type(res_name, &Res_Type);
            polyP->acAm[resn - offset].acAmType = Res_Type;
            gsl_vector_set(polyP->acAm[resn - offset].ca, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].ca, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].ca, 2, z);
            polyP->acAm[resn - offset].caFlag = 1;
        } else if (!strcmp(att, "N")) {
            gsl_vector_set(polyP->acAm[resn - offset].n_, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].n_, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].n_, 2, z);
            polyP->acAm[resn - offset].n_Flag = 1;
        } else if (!strcmp(att, "C")) {
            gsl_vector_set(polyP->acAm[resn - offset].c_, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].c_, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].c_, 2, z);
            polyP->acAm[resn - offset].c_Flag = 1;
        } else if (!strcmp(att, "HN") || !strcmp(att, "H")) {
            gsl_vector_set(polyP->acAm[resn - offset].h_, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].h_, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].h_, 2, z);
            polyP->acAm[resn - offset].h_Flag = 1;
        } else if (!strcmp(att, "O")) {
            gsl_vector_set(polyP->acAm[resn - offset].o_, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].o_, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].o_, 2, z);
            polyP->acAm[resn - offset].o_Flag = 1;
        } else if (!strcmp(att, "CB")) {
            gsl_vector_set(polyP->acAm[resn - offset].cb, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].cb, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].cb, 2, z);
            polyP->acAm[resn - offset].cbFlag = 1;
        } else if (!strcmp(att, "HA")) {
            gsl_vector_set(polyP->acAm[resn - offset].ha, 0, x);
            gsl_vector_set(polyP->acAm[resn - offset].ha, 1, y);
            gsl_vector_set(polyP->acAm[resn - offset].ha, 2, z);
            polyP->acAm[resn - offset].haFlag = 1;
        }
    }
    fclose(fileIn);
    fileIn = NULL;

    return polyP;
}

void WritePdbFile(char *fileName, polyPeptide_t * polyP) {
    int             i = 1, l;
    int             windowSize = polyP->windowSize;
    int             windowOffset = polyP->windowOffset;
    char            nameOut[1000], Resname[4];
    FILE           *FileOut;

    sprintf(nameOut, "%s_%d_%d.pdb", fileName, polyP->acAm[windowOffset].resNb, polyP->acAm[windowSize + windowOffset].resNb);
    FileOut = fopen(nameOut, "w");
    fprintf(FileOut, "\n");
    printf("     fileName : %s \n", nameOut);

    for (l = windowOffset; l <= windowOffset + windowSize; l++) {
        AcAmType2Name(polyP->acAm[l].acAmType, Resname);

        if (polyP->acAm[l].n_Flag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " N  ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].n_), 1., 0., "N");
            i++;
        }
        if (polyP->acAm[l].h_Flag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " H  ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].h_), 1., 0., "H");
            i++;
        }
        if (polyP->acAm[l].caFlag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " CA ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].ca), 1., 0., "C");
            i++;
        }
        if (polyP->acAm[l].haFlag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " HA ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].ha), 1., 0., "H");
            i++;
        }
        if (polyP->acAm[l].cbFlag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " CB ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].cb), 1., 0., "C");
            i++;
        }
        if (polyP->acAm[l].c_Flag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " C  ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].c_), 1., 0., "C");
            i++;
        }
        if (polyP->acAm[l].o_Flag) {
            fprintf(FileOut, STRUCT_FORMAT, "ATOM  ", i, " O  ", "", Resname, "", polyP->acAm[l].resNb, "", XYZ(polyP->acAm[l].o_), 1., 0., "O");
            i++;
        }
    }

    fclose(FileOut);
}

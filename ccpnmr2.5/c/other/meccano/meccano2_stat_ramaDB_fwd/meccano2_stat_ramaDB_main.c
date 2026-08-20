#include "../inc/meccano.h"

int main(int argc, char *argv[]) {
    char            mediaFileName[1000];
    char            csiFileName[1000];
    char            talosFileName[1000];
    char            sequenceFileName[1000];
    char            hbscFileName[1000];
    char            ramaDBFileName[1000];
    char            outputFileName[1000];

    int             firstPPlaneFrag;
    int             lastPPlaneFrag;
    int             ppNbMin;
    int             minValueBest;
    int             maxValueBest;

    int arg = 1;

    if (argv[arg]) {
        strcpy(outputFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(mediaFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(csiFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(talosFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(hbscFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(sequenceFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        strcpy(ramaDBFileName, argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        firstPPlaneFrag = atoi(argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        lastPPlaneFrag = atoi(argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        ppNbMin = atoi(argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        minValueBest = atoi(argv[arg]);
        arg++;
    }
    else
        exit(0);

    if (argv[arg]) {
        maxValueBest = atoi(argv[arg]);
        arg++;
    }
    else
        exit(0);

    RunMeccanoFromFiles(outputFileName, mediaFileName, csiFileName, talosFileName,
        hbscFileName, sequenceFileName, ramaDBFileName, firstPPlaneFrag,
        lastPPlaneFrag, ppNbMin, minValueBest, maxValueBest);

    /*****************************************/
    return 0;
}

#include "phi4.h"

static double M, absM, M2, M4;
static int counter_observables;
static int counter_print = 0;

void initialize_observables () {
    M = 0.; absM = 0.; M2 = 0.; M4 = 0.;
    counter_observables = 0;
}

void update_observables () {
    double auxM = 0.;
    int i;
    for (i = 0; i < V; i++) auxM += phi[i];
    M += auxM;
    absM += fabs(auxM);
    M2 += auxM*auxM;
    M4 += auxM*auxM*auxM*auxM;
    counter_observables ++;
}

void compute_observables () {
    M = M/counter_observables;
    absM = absM/counter_observables;
    M2 = M2/counter_observables;
    M4 = M4/counter_observables;
}

void print_observables (FILE* outFile) {
    counter_print ++;
    fprintf(outFile,"%d \t %f \t %f \t %f \t %f \n", counter_print, M, absM, M2, M4);
}

void get_observables (int* counter_observables_ext, double* mean_m, double* mean_abs_m, double* mean_m2, double* mean_m4) {
    *counter_observables_ext = counter_observables;
    *mean_m = M;
    *mean_abs_m = absM;
    *mean_m2 = M2;
    *mean_m4 = M4;
}


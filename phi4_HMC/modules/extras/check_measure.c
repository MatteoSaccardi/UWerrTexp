#include "phi4.h"

static double DELTAH,ABSDELTAH,EXPDELTAH;
static int COUNTER;
static int COUNTER_PRINT = 0;

void check_initialize_observables () {
    DELTAH = 0.;
    ABSDELTAH = 0.;
    EXPDELTAH = 0.;
    COUNTER = 0;
}

void check_update_observables (double value) {
    DELTAH += value;
    ABSDELTAH += fabs(value);
    EXPDELTAH += exp(-value);
    COUNTER ++;
}

void check_compute_observables () {
    DELTAH /= COUNTER;
    ABSDELTAH /= COUNTER;
    EXPDELTAH /= COUNTER;
}

void check_print_observables (FILE* outFile) {
    COUNTER_PRINT ++;
    fprintf(outFile,"%d \t %f \t %f \t %f \n", COUNTER_PRINT, DELTAH, ABSDELTAH, EXPDELTAH);
}


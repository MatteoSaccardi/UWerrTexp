#ifndef MEASURE_H
#define MEASURE_H

#include "phi4.h"

void initialize_observables ();
void update_observables ();
void compute_observables ();
void print_observables (FILE* outFile);
void get_observables (int* counter_observables_ext, double* mean_m, double* mean_abs_m, double* mean_m2, double* mean_m4);

#endif

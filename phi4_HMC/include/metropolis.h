#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "phi4.h"
#include "measure.h"


/* METROPOLIS_C */
extern double action (act_params_t *apars);
void accept_reject (int *nAccepted, int *nRejected, float kappa, float lambda, float delta, int j);
void check_deltaS (double dS, double r, float kappa, float lambda, float delta, int j);
extern double metropolis (act_params_t *apars, metro_params_t *mpars);
extern double metropolis_with_file (act_params_t *apars, metro_params_t *mpars, FILE* outFile);
extern double metropolis_with_histo (act_params_t *apars, metro_params_t *mpars);
extern double metropolis_prova_binning_naccu (act_params_t *apars, metro_params_t *mpars, FILE* outFile, int naccu, double* mag, double* errMag);




#endif

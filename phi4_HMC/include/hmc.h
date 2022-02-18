#ifndef HMC_H
#define HMC_H

#include "phi4.h"
#include "measure.h"
#include "check_measure.h"
#include "metropolis.h"
#include "istogrammi.h"


/* 3.1 Momentum refreshment */
/* static void initialize_momenta (); */

/* 3.2 Hamiltonian */
/* static double hamiltonian (act_params_t *apars); */

/* 3.3 Molecular dynamics integration */
extern double hmc_with_file (act_params_t *apars, hmc_params_t *hpars, FILE* outFile);
extern double hmc (act_params_t *apars, hmc_params_t *hpars);
/* static void move_phi (double eps);
static double computeForce (act_params_t *apars);
static void move_mom (act_params_t *apars, double eps); */
extern void update_hmc (act_params_t *apars, hmc_params_t *hpars);


/* checks */
extern void check_momenta ();
extern void check_reversibility (act_params_t *apars, hmc_params_t *hpars, double* err_phi, double* err_mom, double* err_H);
extern void hamiltonian_conservation (act_params_t *apars, hmc_params_t *hpars);
extern void checks_analyze (act_params_t *apars, hmc_params_t *hpars);
extern void inexact_vs_exact (act_params_t *apars, hmc_params_t *hpars);

/* prove pre-simulazioni */
extern void PAerfc (act_params_t *apars, hmc_params_t *hpars);


#endif

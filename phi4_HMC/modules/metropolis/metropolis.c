
/*
 *  File metropolic.c
 *  Metropolis algorithm for the phi4 theory
 *
 *  double action(act_params_t *apars)
 *      This routine computes the action S[phi] for the global field phi in
 *      lattice.h and the parameters kappa and lambda from apars.
 *      S = Sum_x [ -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2 ]
 */

#include "metropolis.h"
#include "istogrammi.h"

double action (act_params_t *apars)
{
    int i,mu;
    double phin,S,phi2;
    double kappa = apars->kappa;
    double lambda = apars->lambda;

    S = 0;

    /* loop over all sites */
    for (i = 0; i < V; i++)
    {
        /*sum over neighbors in positive direction*/
        phin = 0;
        for (mu = 0; mu < D; mu++) phin += phi[hop[i][mu]];
        
        phi2 = phi[i]*phi[i];
        S += -2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
    }

    return S;
}

/********************************************************************** 
 * metropolis()
 * Does ntherm+nsweep sweeps over the whole lattice of the Metropolis
 * algorithm. Measurement after each sweep>=ntherm, the averaged measured 
 * values are printed out in fixed intervals, controlled by naccu.
 **********************************************************************/

void check_deltaS (double dS, double r, float kappa, float lambda, float delta, int j)
{
    int mu,i;
    double S = 0., S1 = 0.,phin,phi2;
    for (i = 0; i < V; i++) {
       phin = 0;
       for (mu = 0; mu < D; mu++) phin += phi[hop[i][mu]];
       phi2 = phi[i]*phi[i];
       S += -2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
    }
    phi[j] += delta*(r-0.5); /* new configuration */
    for (i = 0; i < V; i++) {
       phin = 0;
       for (mu = 0; mu < D; mu++) phin += phi[hop[i][mu]];
       phi2 = phi[i]*phi[i];
       S1 += -2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
    }
    phi[j] -= delta*(r-0.5); /* back to starting configuration */
    if ( fabs(S1-S-dS) > 0.00001 ) {
        printf("check_deltaS: il valore teorico %f non coincide con quello calcolato %f \n Premere enter per continuare...",dS,S1-S);
        getchar();
    }
}

void accept_reject (int *nAccepted, int *nRejected, float kappa, float lambda, float delta, int j)
{
    int mu;
    double sumphi = 0.0, phi1 = phi[j], phi2 = phi1*phi1, deltaS, probAcc;
    double r [2];
    double c;
    ranlxd(r,2);
    c = delta*(r[0]-0.5);
    for (mu = 0; mu < 2*D; mu++) sumphi += phi[hop[j][mu]];
    deltaS = -2.0*kappa*c*sumphi + c*(c+2.0*phi1) + lambda*(c*c*(c+2.0*phi1)*(c+2.0*phi1)+2*(phi2-1.0)*c*(c+2*phi1));
    /* check */
    /* check_deltaS(deltaS,r[0],kappa,lambda,delta,j); */
    if (deltaS > 0) probAcc = exp(-1.0*deltaS);
    else probAcc = 1.0;
    if (r[1] < probAcc) {
        phi[j] += delta*(r[0]-0.5);
        *nAccepted += 1;
    }
    else *nRejected += 1;
}

extern double metropolis_with_file (act_params_t *apars, metro_params_t *mpars, FILE* outFile)
{
    /* define parameters */
    int nAccepted = 0;
    int nRejected = 0;

    float kappa = apars->kappa;
    float lambda = apars->lambda;
    int ntherm = mpars->ntherm;
    int nsweep = mpars->nsweep;
    float delta = mpars->delta;
    int naccu = mpars->naccu;
    
    int i, j;
        
    for (i = 0; i < ntherm; i++) {
        for (j = 0; j < V; j++) accept_reject(&nAccepted,&nRejected,kappa,lambda,delta,j);
    }
    /* Ha termalizzato */
    nAccepted = 0;
    nRejected = 0;
    initialize_observables();

    for (i = 0; i < nsweep; i++) {
        for (j = 0; j < V; j++) accept_reject(&nAccepted,&nRejected,kappa,lambda,delta,j);
        update_observables();
        if ((i+1)%naccu == 0) {
            compute_observables();
            print_observables(outFile);
            initialize_observables();
        }
    }

    return (double)nAccepted/(nAccepted+nRejected);
}

double metropolis (act_params_t *apars, metro_params_t *mpars)
{
    FILE* outFile = fopen ("osservabiliMetropolis_default.txt","w");
    double result = metropolis_with_file(apars,mpars,outFile);
    fclose(outFile);
    return result;
}


double metropolis_with_histo (act_params_t *apars, metro_params_t *mpars)
{
    /* define parameters */
    int nAccepted = 0;
    int nRejected = 0;

    float kappa = apars->kappa;
    float lambda = apars->lambda;
    int ntherm = mpars->ntherm;
    int nsweep = mpars->nsweep;
    float delta = mpars->delta;
    int naccu = mpars->naccu;
    
    int i, j;
    
    FILE* outFile = fopen ("defaultEvoluzioneMetropolis.txt","w");
    
    int nValues = nsweep/naccu;
    double* values = (double*) malloc(nValues*sizeof(double));
    
    int counter;
    double m,abs_m,m2,m4;
        
    for (i = 0; i < ntherm; i++) {
        for (j = 0; j < V; j++) accept_reject(&nAccepted,&nRejected,kappa,lambda,delta,j);
    }
    /* Ha termalizzato */
    nAccepted = 0;
    nRejected = 0;
    initialize_observables();

    for (i = 0; i < nsweep; i++) {
        for (j = 0; j < V; j++) accept_reject(&nAccepted,&nRejected,kappa,lambda,delta,j);
        update_observables();
        if ((i+1)%naccu == 0) {
            compute_observables();
            print_observables(outFile);
            get_observables(&counter,&m,&abs_m,&m2,&m4);
            values[(i+1)/naccu-1] = m/V;
            initialize_observables();
        }
    }

    fclose(outFile);
    printf("Su gnuplot, scrivere e.g.\nplot [0:200] 'evoluzioneOsservabili.txt' using 1:2 with linespoints \n");
    create_histogram(100,values,nValues,"istogrammaMagnetizzazione.txt");
    /* create_histogram_in_interval(100,values,nValues,"istogrammaMagnetizzazione.txt",0.5,1.0); */
    /* return (double)nAccepted/(nsweep+ntherm) */
    /* cambia pochissimo */
    return (double)nAccepted/(nAccepted+nRejected);
}

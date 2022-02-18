#include "hmc.h"

static double phi_old[V];
static int NACC = 0;
static int NREJ = 0;
static int CHECKMODE = 0;
static double DELTAH_CHECK = 0.;
static double EXPDELTAH_CHECK = 0.;
static int COUNTER_CHECK = 0;
static FILE* checkFile;



static void initialize_momenta () {
    double r[2];
    int i;
    for (i = 0; i < V/2; i++) {
        ranlxd(r,2);
        mom[2*i] = sqrt(-2.*log(1.-r[0]))*cos(2*M_PI*(1.-r[1]));
        mom[2*i+1] = sqrt(-2.*log(1.-r[0]))*sin(2*M_PI*(1.-r[1]));
    }
    if ((V%2) != 0) {
        ranlxd(r,2);
        mom[V] = sqrt(-2.*log(1.-r[0]))*cos(2*M_PI*(1.-r[1]));
    }
}

static double hamiltonian (act_params_t *apars) {
    double pi2 = 0.,S;
    int i;
    S = action(apars);
    for (i = 0; i < V; i++) pi2 += mom[i]*mom[i];
    return (S+pi2/2.);
}

static void move_phi (double eps) {
    int i;
    for (i = 0; i < V; i++) phi[i] += eps*mom[i];
}

static double computeForce (act_params_t *apars, int i) {
    double sumphi = 0.;
    float kappa = apars->kappa;
    float lambda = apars->lambda;
    int mu;
    for (mu = 0; mu < 2*D; mu++) sumphi += phi[hop[i][mu]];
    return (-2.*kappa*sumphi+2.*phi[i]+4.*lambda*phi[i]*(phi[i]*phi[i]-1.));
}

static void move_mom (act_params_t *apars, double eps) {
    int i;
    for (i = 0; i < V; i++) mom[i] -= eps*computeForce(apars,i);
}

void update_hmc (act_params_t *apars, hmc_params_t *hpars) {
    int i;
    double eps = (hpars->tlength)/(hpars->nstep);

    /* evoluzione di tlength */
    move_mom(apars,eps/2);
    for (i = 0; i < hpars->nstep-1; i++) {
        move_phi(eps);
        move_mom(apars,eps);
    }
    move_phi(eps);
    move_mom(apars,eps/2);
    
    
}

static void accept_reject_hmc (act_params_t *apars, hmc_params_t *hpars) {
    int i;
    double r[1];
    double deltaH;
    initialize_momenta();
    deltaH = hamiltonian(apars);
    for (i = 0; i < V; i++) phi_old[i] = phi[i];
    update_hmc(apars,hpars);
    deltaH = hamiltonian(apars)-deltaH;
    if (CHECKMODE) {
        check_update_observables(deltaH);
        DELTAH_CHECK += fabs(deltaH);
        EXPDELTAH_CHECK += exp(-deltaH);
        COUNTER_CHECK ++;
    }

    if (deltaH > 0.) {
        ranlxd(r,1);
        if (r[0] > exp(-deltaH)) { /* rejected */
            for (i = 0; i < V; i++) phi[i] = phi_old[i];
            NREJ += 1;
        }
        else NACC += 1;
    }
    else NACC += 1;
}

double hmc_with_file (act_params_t *apars, hmc_params_t *hpars, FILE* outFile) {
    int i;
    
    /* ntherm steps di termalizzazione */
    NACC = 0;
    NREJ = 0;
    if (CHECKMODE) {
        COUNTER_CHECK = 0;
        DELTAH_CHECK = 0.;
        EXPDELTAH_CHECK = 0.;
    }
    for (i = 0; i < hpars->ntherm; i++) accept_reject_hmc(apars,hpars);
    
    /* una volta termalizzato, calcolo osservabili */
    initialize_observables();
    NACC = 0;
    NREJ = 0;
    if (CHECKMODE) {
        COUNTER_CHECK = 0;
        DELTAH_CHECK = 0.;
        EXPDELTAH_CHECK = 0.;
        check_initialize_observables();
    }
    for (i = 0; i < hpars->ntraj; i++) {
        accept_reject_hmc(apars,hpars);
        update_observables();
        if ( ((i+1)%hpars->naccu) == 0) {
            compute_observables();
            print_observables(outFile);
            initialize_observables();
            if (CHECKMODE) {
                check_compute_observables();
                check_print_observables(checkFile);
                check_initialize_observables();
            }
        }
    }
    if (CHECKMODE) {
        DELTAH_CHECK /= COUNTER_CHECK;
        EXPDELTAH_CHECK /= COUNTER_CHECK;
    }
    
    return (double)NACC/(NACC+NREJ);
}

double hmc (act_params_t *apars, hmc_params_t *hpars) {
    FILE* outFile = fopen("defaultEvoluzioneHMC.txt","w");
    double accRate = hmc_with_file(apars,hpars,outFile);
    fclose(outFile);
    return accRate;
}


/* checks */

extern void check_momenta () {
    /* CHECKMODE = 1; non serve qui*/
    int N = 1000,i,j,nBins;
    double* momenta = (double*) malloc(N*V*sizeof(double));
    for (i = 0; i < N; i++) {
        initialize_momenta();
        for (j = 0; j < V; j++) momenta[V*i+j] = mom[j];
    }
    nBins = 40;
    create_histogram(nBins, momenta, N*V, "istogrammaMomenti.txt");
}

/* check 1 */
extern void check_reversibility (act_params_t *apars, hmc_params_t *hpars, double* err_phi, double* err_mom, double* err_H) {
    /* CHECKMODE = 1; non serve qui*/
    int i;
    double* mom_old = (double*) malloc(V*sizeof(double));
    initialize_momenta();
    *err_H = hamiltonian(apars);
    for (i = 0; i < V; i++) {
        phi_old[i] = phi[i];
        mom_old[i] = mom[i];
    }
    update_hmc(apars,hpars);
    for (i = 0; i < V; i++) mom[i] = -1.*mom[i];
    update_hmc(apars,hpars);
    for (i = 0; i < V; i++) mom[i] = -1.*mom[i];
    *err_phi = 0.;
    *err_mom = 0.;
    for (i = 0; i < V; i++) {
        *err_phi += fabs(phi_old[i]-phi[i]);
        *err_mom += fabs(mom_old[i]-mom[i]);
    }
    *err_phi /= V;
    *err_mom /= V;
    *err_H = hamiltonian(apars)-(*err_H);
    printf("\no-------------------------------------------------------o\n");
    printf("Check reversibility with tau0 = %f, nstep = %d\n",hpars->tlength,hpars->nstep);
    printf("err_phi = %20.16e, err_mom = %20.16e, err_H = %20.16e\n",*err_phi,*err_mom,*err_H);
    printf("o-------------------------------------------------------o\n\n");

}

/* check 2 */
extern void hamiltonian_conservation (act_params_t *apars, hmc_params_t *hpars) {
    /* stima SENZA analisi degli errori del check 2 */
    FILE* outFile = fopen("hamiltonianConservation.txt","w");
    int N = 5,i;
    double eps,accRate;
    
    CHECKMODE = 1;
    checkFile = fopen("checkFile_hamiltonian.txt","w");
    printf("\no-----------------------o\nhamiltonian_conservation\ntlength = %f\nkappa = %f\nlambda = %f\no-----------------------o\n\n",hpars->tlength,apars->kappa,apars->lambda);
    hpars->ntherm = 100000;
    hpars->ntraj = 100000;
    
    for (i = 0; i < 10; i++) {
        eps = hpars->tlength/N;
        hpars->nstep = N;
        ranlxd(phi,V);
        accRate = hmc(apars,hpars);
        printf("eps = %f, ACCRATE = %f\n",eps,accRate);
        fprintf(outFile,"%f \t %f \t %f\n",eps*eps,DELTAH_CHECK,EXPDELTAH_CHECK);
        N = N+3;
    }
    printf("Su gnuplot, scrivi:\nplot 'hamiltonianConservation.txt' u 1:2\n");
    fclose(outFile);
    fclose(checkFile);
}

/* check 2 e check 4 per analisi errori */
extern void checks_analyze (act_params_t *apars, hmc_params_t *hpars) {

    int N = 5,i;
    double eps,accRate;
    
    CHECKMODE = 1;
    printf("\no-----------------------o\nchecks_analyze\ntlength = %f\nkappa = %f\nlambda = %f\no-----------------------o\n\n",hpars->tlength,apars->kappa,apars->lambda);
    hpars->ntherm = 10000;
    hpars->ntraj = 100000; /* 10000000 per avere più statistica */
    hpars->naccu = 1000;
    checkFile = fopen("checks_analyze.txt","w");
    fprintf(checkFile,"%% N = %d\n",hpars->ntraj/hpars->naccu);
    for (i = 0; i < 10; i++) {
        eps = hpars->tlength/N;
        fprintf(checkFile,"%% eps = %f\n",eps);
        hpars->nstep = N;
        ranlxd(phi,V);
        accRate = hmc(apars,hpars);
        printf("eps = %f, ACCRATE = %f\n",eps,accRate);
        N = N+2;
    }
    
    fclose(checkFile);
}

/* check 3: inexact vs exact algorithm */
extern void inexact_vs_exact (act_params_t *apars, hmc_params_t *hpars) {
    FILE* outFile = fopen("inexact_vs_exact.txt","w");
    int i,j;
    double eps,accRate;
    apars->kappa = 0.18169;
    apars->lambda = 1.3282;
    hpars->tlength = 1.0000;
    hpars->ntherm = 10000;
    hpars->ntraj = 100000; /* 10000000 per avere più statistica */
    hpars->naccu = 1000;
    /* dall'analisi della funzione di correlazione, vediamo che l'autocorrelazione va quasi subito a zero già con 5 step */
    hpars->nstep = 5;
    printf("\no-----------------------o\ninexact_vs_exact\ntlength = %f\nkappa = %f\nlambda = %f\no-----------------------o\n\n",hpars->tlength,apars->kappa,apars->lambda);
    fprintf(outFile,"%% N = %d\n\n",hpars->ntraj/hpars->naccu);
    for (i = 0; i < 10; i++) {
        /* con accept_reject */
        ranlxd(phi,V);
        eps = (hpars->tlength)/(hpars->nstep);
        fprintf(outFile,"%%INIZIO eps = %f\n",eps);
        /* con accept_reject */
        fprintf(outFile,"%%Accept & Reject\n");
        accRate = hmc_with_file(apars,hpars,outFile);
        printf("eps = %f, ACCRATE = %f\n",eps,accRate);
        fprintf(outFile,"%%ACCRATE = %f\n\n",accRate);
        /* senza accept_reject */
        fprintf(outFile,"%%NO Accept & Reject\n");
        ranlxd(phi,V);
        for (j = 0; j < hpars->ntherm; j++) {
            initialize_momenta();
            update_hmc(apars,hpars);
        }
        initialize_momenta();
        for (j = 0; j < hpars->ntraj; j++) {
            initialize_momenta();
            update_hmc(apars,hpars);
            update_observables();
            if ((j+1)%hpars->naccu == 0) {
                compute_observables();
                print_observables(outFile);
                initialize_observables();
            }
        }
        fprintf(outFile,"%%FINE eps = %f\n\n",eps);
        hpars->nstep += 3;
    }
}

/*  */

extern void PAerfc (act_params_t *apars, hmc_params_t *hpars) {
    int i;
    double accRate;
    FILE* PACCfile = fopen("PACCDATAPAerfc.txt","w");
    
    fprintf(PACCfile,"%%File contenente le probabilità di accettanza\n\n");
    CHECKMODE = 1;
    checkFile = fopen("PAerfc.txt","w");
    fprintf(checkFile,"%%SIMULAZIONI PER PROVARE LA RELAZIONE PA = erfc[0.5*sqrt(<dH>)]\n");
    fprintf(checkFile,"%%VOLUME = %d\n",V);
    hpars->nstep = 5;
    apars->kappa = 0.18169;
    apars->lambda = 1.145;
    hpars->tlength = 1.0000;
    hpars->ntherm = 10000;
    hpars->ntraj = 100000;
    hpars->naccu = 1000;
    fprintf(checkFile,"%%MISURA <dH> OGNI N = %d\n",hpars->ntraj/hpars->naccu);
    for (i = 0; i < 10; i++) {
        fprintf(checkFile,"%%deltaTau = %f\n\n",hpars->tlength/hpars->nstep);
        accRate = hmc(apars,hpars);
        fprintf(checkFile,"%%ACCRATE = %f\n\n",accRate);
        fprintf(PACCfile,"%f\n",accRate);
        hpars->nstep = hpars->nstep+2;
    }
    
    fclose(checkFile);
    fclose(PACCfile);
}

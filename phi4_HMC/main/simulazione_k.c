
/* 
 *   File phi4.c
 *
 *   Contains the main program and a few other routines from which the
 *   code to simulate the phi**4 theory can be built. Routines for reading 
 *   in the main parameters of the action and the algorithm are provided.
 *
 *   static int get_val(FILE* fp, char *str, char* fmt,  void* val)
 *      Routine which reads one line from the input file.
 *      Format of the lines is <keyword> <value>.
 *      Checks if the keyword in string str matches,
 *      then gets the value according to the format in fmt
 *
 *   static int read_input(char *input)
 *      Parses the input file (format as specified in get_val)
 *      and prints the parameters onto the screen. Currently
 *      it reads the basic values for the action and also for the 
 *      future Metropolis and the seed of the random number generator.
 */ 

#define CONTROL
#include "phi4.h"
#include "metropolis.h"
#include "hmc.h"

/*  
 *  data structures to store all the parameters of the algorithm,
 *  and action defined in phi4.h
 *  seed for initializing the random numbers
 */

static metro_params_t metro_params;
static act_params_t act_params;
static int seed;
static hmc_params_t hmc_params;
/* static const double DESIREDPACC = 0.9000; */


static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
    char c[128];
    
    if (1 != fscanf(fp,"%s",c)) {
        fprintf(stderr,"Error reading input file at %s\n",str);
        exit(1);
    }
    
    if (strcmp(str,c) != 0) {
        fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
        exit(1);
    }
    
    if (1 != fscanf(fp,fmt,val))
    {
        fprintf(stderr,"Error reading input file at %s\n",str);
        fprintf(stderr,"Cannot read value format %s\n",fmt);
        exit(1);
    }
    
    return 0;
}


static int read_input (char *input, FILE* outFile)
{
    FILE* fp;
    
    fp = fopen(input,"r");
    
    if (fp == NULL)
    {
        fprintf(stderr, "Cannot open input file %s \n",input);
        exit(1);

    }
    
    get_val(fp, "kappa",       "%lf",&act_params.kappa    );
    get_val(fp, "lambda",      "%lf",&act_params.lambda   );
    get_val(fp, "ntherm",      "%i" ,&metro_params.ntherm );
    get_val(fp, "nsweep",      "%i" ,&metro_params.nsweep );
    get_val(fp, "delta",       "%lf",&metro_params.delta  );
    get_val(fp, "seed",        "%i" ,&seed                );
    get_val(fp, "naccu",       "%i" ,&metro_params.naccu  );
    get_val(fp, "tlength",     "%lf",&hmc_params.tlength  );
    get_val(fp, "nstep",       "%i" ,&hmc_params.nstep    );

    
    hmc_params.ntherm = metro_params.ntherm;
    hmc_params.ntraj = metro_params.nsweep;
    hmc_params.naccu = metro_params.naccu;
    
    fprintf(outFile,"%% PARAMETERS       \n"                      );
    fprintf(outFile,"%% L              %i\n", L                   );
    fprintf(outFile,"%% DIM            %i\n", D                   );
    fprintf(outFile,"%% kappa          %f\n", act_params.kappa    );
    fprintf(outFile,"%% lambda         %f\n", act_params.lambda   );
    fprintf(outFile,"%% ntherm         %i\n", metro_params.ntherm );
    fprintf(outFile,"%% nsweep         %i\n", metro_params.nsweep );
    fprintf(outFile,"%% delta          %f\n", metro_params.delta  );
    fprintf(outFile,"%% naccu          %i\n", metro_params.naccu  );
    fprintf(outFile,"%% tlength        %f\n", hmc_params.tlength  );
    fprintf(outFile,"%% nstep          %i\n", hmc_params.nstep    );
    fprintf(outFile,"%% END PARAMETERS   \n"                      );
    
    printf("PARAMETERS       \n"                      );
    printf("L              %i\n", L                   );
    printf("DIM            %i\n", D                   );
    printf("kappa          %f\n", act_params.kappa    );
    printf("lambda         %f\n", act_params.lambda   );
    printf("ntherm         %i\n", metro_params.ntherm );
    printf("nsweep         %i\n", metro_params.nsweep );
    printf("delta          %f\n", metro_params.delta  );
    printf("naccu          %i\n", metro_params.naccu  );
    printf("tlength        %f\n", hmc_params.tlength  );
    printf("nstep          %i\n", hmc_params.nstep    );
    printf("END PARAMETERS   \n"                      );
    
    return 0;
}

void printMATLABCommentLine (FILE* outFile) {
    fprintf(outFile,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
}

int main (int argc, char* argv[])
{
    double accRate,deltaTau;
    FILE* outFile = fopen("parametriInfile.txt","w");
    char nameFile [150];

    if (argc != 2)
    {
       fprintf(stderr,"Number of arguments not correct\n");
       fprintf(stderr,"Usage: %s <infile> \n",argv[0]);
       exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1],outFile);
    fclose(outFile);

    /* initialize random number generator */
    rlxd_init(1,seed);

    /* initialize the nearest neighbor field */
    hopping(hop);

    /* initialize phi field */
    ranlxd(phi,V);
    
    /* D, L, V sono fissati */
    
    /* Dato che vogliamo una probabilità di accettanza circa costante e pari a DESIREDPACC, e dato che abbiamo verificato che il prodotto V*(deltaTau)^4 deve restare circa costante affinché PACC lo sia altrettanto, in base ai valori fissati di D, L, V calcoliamo deltaTau, da cui il numero di step (tlenght = 1 fissato). */
    /* Nota: abbiamo trovato il valore DESIREDPACC per nStep = 8, V = 64. */
    deltaTau = (1./8)*pow(64./V,1./4);
    hmc_params.nstep = hmc_params.tlength/deltaTau; /* è un intero, approssima il valore esatto */
    
    /* Le simulazioni che facciamo sono a lambda = 1.145 fissato. */
    if (act_params.lambda != 1.145) {
        printf("Le nostre simulazioni sono a lambda = 1.145 fissato, ma da infile sembra che lambda sia diverso. Lo riportiamo al valore desiderato.\n");
        act_params.lambda = 1.145;
    }
    printf("\nD = %d, L = %d, V = %d\n\n",D,L,V);
    
    /* Simulazione */
    sprintf(nameFile,"simulazione_k%.0f.txt",act_params.kappa*1e6);
    outFile = fopen(nameFile,"a");
    printMATLABCommentLine(outFile);
    fprintf(outFile,"%%SIMULAZIONE D = %d, L = %d, V = %d\n",D,L,V);
    fprintf(outFile,"%%NUMERO DI MISURE PER OGNI KAPPA: T = %d\n",hmc_params.ntraj/hmc_params.naccu);
    printMATLABCommentLine(outFile);
    printf("Kappa = %f\n",act_params.kappa);
    accRate = hmc_with_file(&act_params,&hmc_params,outFile);
    printf("ACCRATE = %f\n\n",accRate);
    fprintf(outFile,"%%ACCRATE = %f\n\n",accRate);
    printf("Cambia il valore di L in ../include/lattice.h, fai make e quindi esegui.\n\n");
    
    return 0;
}


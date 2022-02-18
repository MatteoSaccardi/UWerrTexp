
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

    return 0;
}


int main (int argc, char* argv[])
{
    FILE* outFile = fopen("datiInfile.txt","w");
    int i,lato;
    double deltaTau,accRate;

    if (argc != 2)
    {
       fprintf(stderr,"Number of arguments not correct\n");
       fprintf(stderr,"Usage: %s <infile> \n",argv[0]);
       exit(1);
    }

    /* get the parameters from the input file */
    read_input(argv[1],outFile);
    fclose(outFile);
    
    /* here we want to use these parameters */
    act_params.kappa = 0.18169;
    act_params.lambda = 1.145;
    hmc_params.tlength = 1.0000;
    hmc_params.ntherm = 10000;
    hmc_params.ntraj = 100000;
    hmc_params.naccu = 1000;
    /* hmc_params.nstep è fissato da fuori, varia in modo tale da mantenere PACC circa costante; i valori da usare in funzione di L sono stampati (vedi sotto) */

    /* initialize random number generator */
    rlxd_init(1,seed);
    /* initialize the nearest neighbor field */
    hopping(hop);
    /* initialize phi field */
    ranlxd(phi,V);
    
    /* Simulazione */
    accRate = hmc(&act_params,&hmc_params);
    outFile = fopen("provePAcostante.txt","a");
    fprintf(outFile,"%d\t%f\t%f\n",V,pow(hmc_params.tlength/hmc_params.nstep,4),accRate);
    /* vogliamo verificare che V*(deltaTau)^4 costante porti a accRate circa costante */
    fclose(outFile);
    
    /* Con deltaTau = 1/8, V = 64 abbiamo trovato PA = 0.868760. Proviamo a mantenerlo costante */
    /* OSS: V è fissato, devo cambiarlo a mano per fare varia simulazioni */
    printf("ProvePAcostante: partendo da L = 4, nstep = 8, deltaTau = 1/8 varia L e deltaTau (nstep) come segue per verificare che PACC rimanga circa costante. Di seguito stampiamo, in ordine, i valori di L, deltaTau, nStep da usare.\n");
    for (i = 0; i < 7; i++) {
        lato = 4+2*i; /* D = 3 */
        deltaTau = (1./8)*pow(64./pow(lato,3),1./4);
        printf("%d\t%f\t%f\n",lato,deltaTau,hmc_params.tlength/deltaTau);
    }
    printf("Ora, hai fatto la prova con L = %d, nstep = %d\n",L,hmc_params.nstep);
    printf("Hai ottenuto PACC = %f\n",accRate);
    printf("Modifica L in lattice.h e nstep in infile. Dopodiché, rifai make ed esegui.\n");
    printf("Ricordati che sovrascrivi il file provePAcostante.txt, quindi se hai finito salvalo altrove o con un altro nome prima di iniziare diverse simulazioni e cancella questo.\n");
    
    
     

    return 0;
}


#ifndef LATTICE_H
#define LATTICE_H

/* dimension of the lattice */
#define D 3
/* spatial extend of the lattice */
#define L 8
/* lattice volume, needs to be adjusted according to number of dimensions */
#define V (L*L*L)

#ifdef CONTROL 
#define EXTERN 
#undef CONTROL
#else
#define EXTERN extern
#endif

EXTERN double phi[V];
EXTERN int    hop[V][2*D];
EXTERN double mom[V];

#endif

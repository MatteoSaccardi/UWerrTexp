#include "phi4.h"

void create_histogram (int nBins, double* values, int nValues, char* fileName) {
    double maxValue = 0., minValue = 0., step = 0.;
    int i, j, maxCounts;
    int* nCounts = (int*) malloc(nBins*sizeof(int));
    FILE* outFile = fopen(fileName,"w");
    
    /* Riempio l'istogramma */
    for (i = 0; i < nBins; i++) nCounts[i] = 0.;
    for (i = 0; i < nValues; i++) {
        if (values[i] > maxValue) maxValue = values[i];
        if (values[i] < minValue) minValue = values[i];
    }
    step = (maxValue-minValue)/nBins;
    maxCounts = 0.;

    for (i = 0; i < nValues; i++) {
        j = (values[i]-minValue)/step;
        nCounts[j]++;
    }
    for (i = 0; i < nBins; i++) {
        if (nCounts[i] > maxCounts) maxCounts = nCounts[i];
    }
    /* Stampo l'istogramma */
    for (i = 0; i < nBins; i++) {
        fprintf(outFile,"%f \t %f \n",minValue+(i+0.5)*step,(double)nCounts[i]/maxCounts);
    }
}

void create_histogram_in_interval (int nBins, double* values, int nValues, char* fileName, double minValue, double maxValue) {
    double step = 0.;
    int i, j, maxCounts;
    int* nCounts = (int*) malloc(nBins*sizeof(int));
    FILE* outFile = fopen(fileName,"w");
    
    /* Riempio l'istogramma */
    for (i = 0; i < nBins; i++) nCounts[i] = 0.;
    step = (maxValue-minValue)/nBins;
    maxCounts = 0.;

    for (i = 0; i < nValues; i++) {
        if (values[i] < minValue || values[i] > maxValue) continue;
        j = (values[i]-minValue)/step;
        nCounts[j]++;
    }
    for (i = 0; i < nBins; i++) {
        if (nCounts[i] > maxCounts) maxCounts = nCounts[i];
    }
    /* Stampo l'istogramma */
    for (i = 0; i < nBins; i++) {
        fprintf(outFile,"%f \t %d \n",minValue+i*step,nCounts[i]);
    }
}

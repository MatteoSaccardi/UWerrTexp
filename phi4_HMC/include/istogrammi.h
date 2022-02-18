#ifndef ISTOGRAMMI_H
#define ISTOGRAMMI_H

#include "phi4.h"

void create_histogram (int nBins, double* values, int nValues, char* fileName);
void create_histogram_in_interval (int nBins, double* values, int nValues, char* fileName, double minValue, double maxValue);

#endif

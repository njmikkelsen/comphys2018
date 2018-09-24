#ifndef JACOBIALGORITHM_H
#define JACOBIALGORITHM_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "time.h"

using namespace std;

double ** getJacobiEigenpairs (const int n, double **A, double tol, int N_max);
void      doGivensSimTrans    (const int n, double **A, double **V, int k, int l);
double    getMaxOffDiag       (const int n, double **A, int *k, int *l);

#endif

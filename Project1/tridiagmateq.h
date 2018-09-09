#ifndef TRIDIAGMATEQ_H
#define TRIDIAGMATEQ_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include "time.h"
#include "armadillo"

using namespace std;
using namespace arma;

double* gen_x (int n, double h);
double* gen_f (int n, double h, double* x);
double* gen_y (int n, double h, double* x);

void save_data (int n, string alg, double dt, double T, double* x, double* exact, double* y);

void arma_tridiagmateq_solver     (int n, double* a, double* b, double* c, double* y, double* f,
                                   clock_t* t0, clock_t* t1);
void general_tridiagmateq_solver  (int n, double* a, double* b, double* c, double* y, double* f,
                                   clock_t* t0, clock_t* t1);
void special_tridiagmateq_solver  (int n, double a, double b, double c, double* y, double* f,
                                   clock_t* t0, clock_t* t1);
void taylored_tridiagmateq_solver (int n, double* y, double* f, clock_t* t0, clock_t* t1);


#endif 

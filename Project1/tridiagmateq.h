#ifndef TRIDIAGMATEQ_H
#define TRIDIAGMATEQ_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

double* gen_x (int n, double h);
double* gen_f (int n, double h, double* x);
double* gen_y (int n, double h, double* x);

int save_data (int n, string algorithm, double* x, double* exact, double* y);
inline bool file_existence (const string& filename);

#endif

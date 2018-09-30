#ifndef MISCLIB_H
#define MISCLIB_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;

void display_matrix (const int n, double **M);
void display_matrix (const int n, double **M, string matrix);
void display_matrix (const int n, double **M,                int precision);
void display_matrix (const int n, double **M, string matrix, int precision);

void display_vector (const int n, double *v);
void display_vector (const int n, double *v, string vector);
void display_vector (const int n, double *v,                int precision);
void display_vector (const int n, double *v, string vector, int precision);

double ** init_matrix                  (const int m, const int n);
double ** init_matrix                  (const int n);
double ** init_random_matrix           (const int m, const int n);
double ** init_random_matrix           (const int n);
double ** init_symmetric_random_matrix (const int m, const int n);
double ** init_symmetric_random_matrix (const int n);
double ** init_identity_matrix         (const int m, const int n);
double ** init_identity_matrix         (const int n);
double ** init_tridiag_Toeplitz_matrix (const int m, const int n, double a, double b, double c);
double ** init_tridiag_Toeplitz_matrix (const int n,               double a, double b, double c);
void      delete_matrix                (const int n, double ** M);

void     scale_matrix        (const int n, double **M, double c);
double * get_matrix_diagonal (const int n, double **M);
double * sort_vector         (const int n, double *v, int *idx);

#endif

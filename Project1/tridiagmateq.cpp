#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <string>
#include "tridiagmateq.h"
#include "armadillo"
#include "time.h"

using namespace std;
using namespace arma;

double* gen_x (int n, double h) {
  int     i;
  double* x = new double[n]; 
  for (i=0; i<n; i++) {x[i] = (i+1)*h;}
  return x;
};

double* gen_f (int n, double h, double* x) {
  int     i;
  double* f = new double[n]; 
  for (i=0; i<n; i++) {f[i] = -100*exp(-10*x[i]);}
  return f;
};

double* gen_y (int n, double h, double* x) {
  int     i;
  double* y = new double[n]; 
  for (i=0; i<n; i++) {y[i] = 1-(1-exp(-10))*x[i]-exp(-10*x[i]);}
  return y;
};

void save_data (int n, string alg, double dt, double T,
                double* x, double* exact, double* y) {
  int i;
  string filename = to_string(n) + "_" + alg + ".dat";
  ofstream outfile("data/"+filename, ios::trunc);
  outfile << "n         = " << n << endl;
  outfile << "algorithm = " << alg << endl;
  outfile << "no. ticks = " << setprecision(16) <<  T << endl;
  outfile << "ticks/sec = " << setprecision(16) << dt << endl;
  outfile << "-----------------------------------------------------" << endl;
  outfile << "        x        |     y_exact     |      y_num" << endl;
  for (i=0; i<n; i++) {
    outfile << " " << left << scientific << setprecision(9) << setw(10)
            << x[i] << " | " << exact[i] << " | " << y[i] << endl;
  }
  cout << "data saved in 'data/" << filename << "'." << endl;
};

void arma_tridiagmateq_solver (int n, double* a, double* b, double* c, double* y, double* f,
                               clock_t* t0, clock_t* t1) {
    int i;
    // armadillo declerations
    mat D(n,n), L(n,n), U(n,n);
    vec f_vec(n), y_vec(n), Uy(n);
    // definitions
    for (i=0; i<n-1; i++) {
      D(i+1,i) = a[i];
      D(i,i)   = b[i];
      D(i,i+1) = c[i];
      f_vec(i) = f[i];
    }
    D(n-1,n-1) = b[n-1];
    f_vec[n-1] = f[n-1];
    // solve matrix equations
    *t0    = clock();
    lu(L,U,D);
    Uy     = solve(L,f_vec);
    y_vec  = solve(U,Uy);
    *t1    = clock();
    // return to y array
    for (i=0; i<n; i++) {y[i] = y_vec(i);}
};

void general_tridiagmateq_solver (int n, double* a, double* b, double* c, double* y, double* f,
                                  clock_t* t0, clock_t* t1) {
  int i;
  double alpha;
  double* c_new = new double[n-1];
  double* f_new = new double[n-1]; // final f_new is directly defined as y[n-1]
  *t0           = clock();
  c_new[0]      = c[0]/b[0];
  f_new[0]      = f[0]/b[0];
  for (i=1; i<n-1; i++) {
    alpha    = b[i] - a[i]*c_new[i-1];
    c_new[i] = c[i]/alpha;
    f_new[i] = (f[i]-a[i]*f_new[i-1])/alpha;
  }
  y[n-1] = (f[n-1]-a[n-1]*f_new[n-2])/(b[n-1]-a[i-1]*c_new[n-2]);
  for (i=n-2; i>=0; i--) {
    y[i] = f_new[i]-c_new[i]*y[i+1];
  }
  *t1 = clock();
};

void special_tridiagmateq_solver  (int n, double a, double b, double c, double* y, double* f,
                                   clock_t* t0, clock_t* t1) {
  int i;
  double rho1, rho2, rho3, alpha;
  double* f_new   = new double[n-1];  // final f_new is directly defined as y[n-1]
  double* c_new   = new double[n-1];
  *t0 = clock();
  rho1 = c/b;
  rho2 = b/a;
  rho3 = c/a;
  f_new[0] = f[0]/b;
  c_new[0] = rho1;
  for (i=1; i<n-1; i++) {
    alpha    = rho2-c_new[i-1];
    c_new[i] = rho3/alpha;
    f_new[i] = (f[i]/a-f_new[i-1])/alpha;
  }
  y[n-1] = (f[n-1]/a-f_new[n-2])/(rho2-c_new[n-2]);
  for (i=n-2; i>=0; i--) {y[i] = f_new[i]-c_new[i]*y[i+1];}
  *t1 = clock();
  delete[] f_new,c_new;
};

void taylored_tridiagmateq_solver (int n, double* y, double* f, clock_t* t0, clock_t* t1) {
  int i;
  double* f_new = new double[n-1];
  double* c_new = new double[n];
  *t0           = clock();
  f_new[0] = -0.5*f[0];
  for (i=1; i<n+1; i++) {c_new[i-1] = (double)i/((double)(i+1));}
  for (i=1; i<n-1; i++) {f_new[i]   = (f_new[i-1]-f[i])*c_new[i];}
  y[n-1] = (f_new[n-2]-f[n-1])*c_new[n-1];
  for (i=n-2; i>=0; i--) {y[i] = f_new[i]+c_new[i]*y[i+1];}
  *t1 = clock();
  delete[] f_new,c_new;
};


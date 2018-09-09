#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include "tridiagmateq.h"
#include "armadillo"
#include "time.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // miscellaneous variable declerations
  int     i, n;
  double  h, hh;
  string  alg;
  clock_t t0, t1;

  // verify command line arguments
  if (argc < 3) {
    cout << "Error: Invalid command line arguments!" << endl;
    cout << "Expected 2 arguments: matrix dim (int) & algorithm (string)" << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }
  n   = atoi(argv[1]);
  alg = argv[2];
  if (n<3) {
    cout << "Error: n must be greater than 2" << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }
  h = 1/((double)n+1);
  
  // generate arrays 
  double* x     = gen_x(n,h);
  double* f     = gen_f(n,h,x);
  double* exact = gen_y(n,h,x);
  double* y     = new double[n];
  double* a     = new double[n-1];
  double* b     = new double[n];
  double* c     = new double[n-1];
  for (i=0; i<n-1; i++) {
    a[i] = 1;
    b[i] = -2;
    c[i] = 1;
  }
  b[n-1] = -2;
  
  // algorithms
  if      (alg == "LU")       {    arma_tridiagmateq_solver(n, a, b, c,  y, f, &t0, &t1);}
  else if (alg == "general")  { general_tridiagmateq_solver(n, a, b, c,  y, f, &t0, &t1);}
  else if (alg == "special")  { special_tridiagmateq_solver(n, 1, -2, 1, y, f, &t0, &t1);}
  else if (alg == "taylored") {taylored_tridiagmateq_solver(n,           y, f, &t0, &t1);}
  
  // Failure to input a valid algorithm name.
  else {
    cout << "Algorithm argument does not match any of the available options:" << endl;
    cout << "LU, general, special, taylored." << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }
  
  // scale solution with h^2
  hh = h*h;
  for (i=0; i<n; i++) {y[i] = hh*y[i];}
  
  // write data to file or print to terminal
  save_data(n, alg, CLOCKS_PER_SEC, t1-t0, x, exact, y);
  
  // free memory
  delete[] x,f,exact,y,a,b,c;

  return 0;
}



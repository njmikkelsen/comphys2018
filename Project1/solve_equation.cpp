#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "tridiagmateq.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // loop variable decleration
  int i;

  // verify command line arguments
  if (argc < 3) {
    cout << "Error: Invalid command line arguments!" << endl;
    cout << "Expected 2 arguments: matrix dim (int) & algorithm (string)" << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }

  // program parameters
  int    n         = atoi(argv[1]);
  string algorithm = argv[2];
  
  // verify n >= 3
  if (n<3) {
    cout << "Error: n must be greater than 2" << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }
  
  // generate arrays
  double  h     = 1/((double)n+1);
  double* x     = gen_x(n,h);
  double* f     = gen_f(n,h,x);
  double* exact = gen_y(n,h,x);
  double* y     = new double[n];
  
  // LU-decomposition
  if (algorithm == "LU") {
    mat D(n,n), L(n,n), U(n,n);
    vec f_vec(n), y_vec(n), Uy(n);
    for (i=0; i<n-1; i++) {
      D(i,i+1) = 1;
      D(i,i)   = -2;
      D(i+1,i) = 1;
      f_vec(i) = f[i];
    }
    D(n-1,n-1) = -2;
    f_vec[n-1] = f[n-1];
    lu(L,U,D);
    Uy     = solve(L,f_vec);
    y_vec  = solve(U,Uy);
    for (i=0; i<n; i++) {
      y[i] = h*h*y_vec(i);
    }
  }
  
  // Generalised tri-diagonal matrix equation algorithm
  else if (algorithm == "general") {
    cout << algorithm << endl;
  }

  // Specialised tri-diagonal matrix equation algorithm
  else if (algorithm == "special") {
    cout << algorithm << endl;
  }

  // Taylored tri-diagonal matrix equation algorithm
  else if (algorithm == "taylored") {
    cout << algorithm << endl;
  }

  // Failure to input a valid algorithm name.
  else {
    cout << "Algorithm argument does not match any of the available options:" << endl;
    cout << "LU, general, special, taylored." << endl;
    cout << "Aborting program." << endl;
    exit(1);
  }
  
  // write data to file
  int success = save_data (n, algorithm, x, exact, y);
  
  // free-up memory
  delete[] x,f,exact,y;
  
  return 0;
}



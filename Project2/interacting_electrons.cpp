#include "jacobi_algorithm.h"
#include "misclib.h"
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // default or command line parameters
  const int    n       = (argc > 1) ? atoi(argv[1]) : 10;    // number of grid points
  const double xi_inf  = (argc > 2) ? atof(argv[2]) : 1e1;   // approximation of infinity
  const double omega   = (argc > 3) ? atof(argv[3]) : 1e0;   // harmonic oscillator frequency
  const int    N_max   = (argc > 3) ? atoi(argv[4]) : 100;   // maximum number of iterations
  const double epsilon = (argc > 4) ? atof(argv[5]) : 1e-8;  // max element error tolerance
  
  // print parameters to file
  ofstream File("InteractingElectrons/"+to_string(n)+"_"+to_string((int)(omega*10.))+".dat", ios::trunc);
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  File << "Numerical Analysis Of Two Coulomb-Interacting Electrons In A Harmonic Oscillator Potential\n\n";
  File << "Program parameters:"      << endl;
  File << "  n        = " << n       << endl;
  File << "  infinity = " << xi_inf  << endl;
  File << "  omega    = " << omega   << endl;
  File << "  N_max    = " << N_max   << endl << scientific;
  File << "  epsilon  = " << epsilon << "\n\n";
  File.flags(f);
  
  // misc
  int i,j;
  int *sort_idx = new int [n];
  vec a(n);   a.fill(-2);
  vec b(n-1); b.fill(1);
  
  // init grid
  double h   = xi_inf/((double)n+1);
  double *xi = new double [n];
  for (i=0; i<n; i++) {xi[i] = (i+1)*h;};
  
  // init D matrix
  double **D = init_tridiag_Toeplitz_matrix(n,1,-2,1);
  for (i=0; i<n; i++) {D[i][i] -= omega*omega*xi[i]*xi[i];};
    
  // redirect cout stream to file
  streambuf * strm_buffer = cout.rdbuf();
  cout.rdbuf(File.rdbuf());
  
  // find eigenpairs
  double **V = getJacobiEigenpairs(n,D,epsilon,N_max);
  cout.rdbuf(strm_buffer);  // redirect back to cout
  double *eigval      = get_matrix_diagonal(n,D);
  double *eigval_sort = sort_vector(n,eigval,sort_idx);
  
  // write eigenvalues and eigenvectors to file
  File << "\nEigenvalues:" << scientific << showpos << endl;
  for (i=0; i<n; i++) {
    File << eigval_sort[i] << endl;
  };
  File << "\nEigenvectors:" << "\n\n";
  for (i=0; i<n; i++) {
    File << "eigenvalue = " << eigval_sort[i] << ":\n";
    for (j=0; j<n; j++) {File << V[j][sort_idx[i]] << endl;};
    File << endl;
  };
  
  // delete D matrix
  delete_matrix(n,D);
  
  return 0;
}

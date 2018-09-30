#include "jacobi_algorithm.h"
#include "misclib.h"
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // default or command line parameters
  const int    n       = (argc > 1) ? atoi(argv[1]) : 100;   // number of grid points
  const int    N_max   = (argc > 2) ? atoi(argv[2]) : 1000;  // maximum number of iterations
  const double epsilon = (argc > 3) ? atof(argv[3]) : 1e-6;  // max element error tolerance
  
  // print parameters to file
  ofstream File("./BucklingBeam/"+to_string(n)+"_"+to_string(-(int)log10(epsilon))+".dat", ios::trunc);
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  File << "Numerical Analysis Of A Buckling Beam" << "\n\n";
  File << "Program parameters:"     << endl;
  File << "  n       = " << n       << endl;
  File << "  N_max   = " << N_max   << endl << scientific;
  File << "  epsilon = " << epsilon << "\n\n";
  File.flags(f);
  
  // misc 
  int i,j;
  double exact;
  int *sort_idx = new int [n];
  vec a(n);   a.fill(-2);
  vec b(n-1); b.fill(1);
  
  // init D matrix
  double **D = init_tridiag_Toeplitz_matrix(n,1,-2,1);
  mat D_arma = diagmat(a) + diagmat(b,-1) + diagmat(b,1);
  
  // redirect cout stream to file
  streambuf * strm_buffer = cout.rdbuf();
  cout.rdbuf(File.rdbuf());
  
  // find eigenpairs
  double **V = getJacobiEigenpairs(n,D,epsilon,N_max);
  cout.rdbuf(strm_buffer);  // redirect back to cout
  double *eigval      = get_matrix_diagonal(n,D);
  double *eigval_sort = sort_vector(n,eigval,sort_idx);
  vec eigval_arma;
  mat eigvec_arma;
  eig_sym(eigval_arma,eigvec_arma,D_arma);
  
  // write eigenvalues to file
  File << "\nEigenvalues:" << endl;
  File << "     Exact    |     Jacobi    |    Armadillo" << scientific << showpos << endl;
  for (i=0; i<n; i++) {
    exact = -2+2*cos((n-i)*M_PI/((double)(n+1)));
    File << exact << " | " << eigval_sort[i] << " | " << eigval_arma[i] << endl;
  };
  
  // write eigenvectors to file
  File << "\nEigenvectors:" << scientific << showpos << endl;
  for (i=0; i<n; i++) {
    File << "\neigenvalue = " << eigval_sort[i] << ":\n";
    File << "    Jacobi    |    Armadillo" << endl;
    for (j=0; j<n; j++) {File << V[j][sort_idx[i]] << " | " << eigvec_arma(i,j) << endl;};
  };
  
  // delete D matrix
  delete_matrix(n,D);
  
  return 0;
}

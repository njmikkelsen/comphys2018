#include "jacobi_algorithm.h"
#include "misclib.h"
#include "armadillo"
using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // default or command line parameters
  const int    n       = (argc > 1) ? atoi(argv[1]) : 10;    // number of grid points
  const int    N_max   = (argc > 2) ? atoi(argv[2]) : 100;   // maximum number of iterations
  const double epsilon = (argc > 3) ? atof(argv[3]) : 1e-8;  // max element error tolerance
  
  // print parameters to command line
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << "Buckling Beam - Testing the Jacobi algorithm" << endl;
  cout << "Program parameters:"     << endl;
  cout << "  n       = " << n       << endl;
  cout << "  N_max   = " << N_max   << endl << scientific;
  cout << "  epsilon = " << epsilon << "\n\n";
  cout.flags(f);
  
  // misc
  int i,j;
  vec a(n);   a.fill(-2);
  vec b(n-1); b.fill(1);
  
  // init D matrix
  double **D = init_tridiag_Toeplitz_matrix(n,1,-2,1);
  mat D_arma = diagmat(a) + diagmat(b,-1) + diagmat(b,1);
  
  // find eigenpairs
  double **V     = getJacobiEigenpairs(n,D,epsilon,N_max);
  double *eigval = get_matrix_diagonal(n,D);
  vec eigval_arma;
  mat eigvec_arma;
  eig_sym(eigval_arma,eigvec_arma,D_arma);
  
  display_matrix(n,D,"Jacobi");
  eigvec_arma.print("Armadillo =");
  
  
  // delete D matrix
  delete_matrix(n,D);
  
  return 0;
}

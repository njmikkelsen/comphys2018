#include "jacobi_algorithm.h"
#include "misclib.h"
#include "armadillo"

using namespace std;
using namespace arma;

/*
command line compilation: (Ubuntu  16.04.5 LTS)
g++ compare_jacobi_armadillo.cpp jacobi_algorithm.cpp misclib.cpp -o compare_jacobi_armadillo.x -std=c++14 -l armadillo

example run:

[terminal]$ ./compare_jacobi_armadillo.x 5 25 1e-12 50
COMPARISON PROGRAM:
JACOBI EIGENPAIR ALGORITHM & ARMADILLO'S STANDARD DIAGONALIZATION ALGORITHM

Program parameters:
  n       = 5
  N_max   = 25
  epsilon = 1.000000e-12
  c       = 5.000000e+01

Jacobi eigenpair algorithm completed:
  Number of iterations = 25
  Max. Abs. off-diag.  = 1.642819e-04
  time spent           = 7.100000e-05 s

Standard Armadillo diagonalization process completed.
  time spent = 1.992570e-04 s

Jacobi eigenvalues =
  -6.177266e+00
  +2.304720e+01
  +3.458100e+00
  -2.783981e+01
  +1.529479e+02
Armadillo eigenvalues =
  -2.7840e+01
  -6.1773e+00
   3.4581e+00
   2.3047e+01
   1.5295e+02
Jacobi eigenvectors = 
  +7.5626e-01    -8.7331e-02    +1.6898e-01    -2.8486e-01    +5.5744e-01  
  -3.2775e-01    +6.1086e-01    -2.1745e-01    -6.2419e-01    +2.8730e-01  
  -4.6547e-01    -9.9409e-02    +7.6965e-01    +7.2050e-02    +4.1942e-01  
  +6.5993e-02    +5.1101e-01    -1.8958e-01    +7.2384e-01    +4.1789e-01  
  -3.1565e-01    -5.9009e-01    -5.4394e-01    +1.0703e-02    +5.0614e-01  
Armadillo eigenvectors =
   0.2849   0.7563  -0.1690  -0.0873   0.5574
   0.6242  -0.3278   0.2174   0.6109   0.2873
  -0.0720  -0.4655  -0.7696  -0.0994   0.4194
  -0.7238   0.0660   0.1896   0.5110   0.4179
  -0.0107  -0.3156   0.5439  -0.5901   0.5061
*/

int main(int argc, char* argv[])
{
  // parameters
  const int    n       = (argc > 1) ? atoi(argv[1]) : 4;;
  const int    N_max   = (argc > 2) ? atoi(argv[2]) : 20;
  const double epsilon = (argc > 3) ? atof(argv[3]) : 1e-8;
  const double c       = (argc > 4) ? atof(argv[4]) : 200;
  
  // misc
  int i,j;
  wall_clock arma_time;
  
  // print parameters
  cout << "COMPARISON PROGRAM:\nJACOBI EIGENPAIR ALGORITHM & ARMADILLO'S STANDARD DIAGONALIZATION ALGORITHM\n\n";
  cout << "Program parameters:"     << endl;
  cout << "  n       = " << n       << endl;
  cout << "  N_max   = " << N_max   << endl << scientific;
  cout << "  epsilon = " << epsilon << endl;
  cout << "  c       = " << c       << "\n\n";
  
  // random square matrix A
  double **A1 = init_symmetric_random_matrix(n);
  scale_matrix(n,A1,c);
  mat A2(n,n);
  for (i=0; i<n; i++) {for (j=0; j<n; j++) {A2(i,j) = A1[i][j];};};
  
  // diagonalize A via Jacobi algorithm
  double ** V1    = getJacobiEigenpairs(n,A1,epsilon,N_max);
  double *lambda1 = get_matrix_diagonal(n,A1);
  cout << "\n";
  
  // diagonalize A via armadillo 
  mat V2(n,n);
  vec lambda2(n);
  arma_time.tic();
  eig_sym(lambda2,V2,A2);
  double t_arma = arma_time.toc();
  cout << "Standard Armadillo diagonalization process completed." << endl;
  cout << "  time spent = " << t_arma << " s\n\n";
  
  
  // print eigenvalues
  cout << "Jacobi eigenvalues =" << showpos << scientific << endl;
  for (i=0; i<n; i++) {cout << "  " << A1[i][i] << endl;};
  lambda2.print("Armadillo eigenvalues =");
  
  // print eigenvectors
  display_matrix(n,V1,"Jacobi eigenvectors");
  V2.print("Armadillo eigenvectors =");
  
  
  // free space
  delete_matrix(n,A1);
  delete_matrix(n,V1);

  return 0;
}

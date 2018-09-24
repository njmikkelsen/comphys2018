#include "jacobi_algorithm.h"
#include "misclib.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main() {
  // misc
  int i,j;

  // parameters
  const int    n       = 4;
  const double epsilon = 1.0e-8;
  const int    N_max   = 20;
  const double c       = 1.0e+2;
  
  // random square matrix A
  double **A1 = init_symmetric_random_matrix(n);
  scale_matrix(n,A1,c);
  mat A2(n,n);
  for (i=0; i<n; i++) {for (j=0; j<n; j++) {A2(i,j) = A1[i][j];};};
  
  // diagonalize A via Jacobi algorithm
  double ** V1    = getJacobiEigenpairs(n,A1,epsilon,N_max);
  double *lambda1 = get_matrix_diagonal(n,A1);
  
  // diagonalize A via armadillo 
  mat V2(n,n);
  vec lambda2(n);
  eig_sym(lambda2,V2,A2);
  
  
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

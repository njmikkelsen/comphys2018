#include "jacobi_algorithm.h"
#include "misclib.h"

double ** getJacobiEigenpairs (const int n, double **A, double tol, int N_max) {
  /*
  Diagonalize and extract the eigenpairs of a real symmetric matrix via the Jacobi
  eigenpair algorithm.
  
  Arguments:
    n     | matrix dimensionality
    A     | matrix to diagonalize
    V     | an identity matrix
    tol   | error tolerance
    N_max | maximum number of iterations
  */

  // misc
  double **V = init_identity_matrix(n);
  int i = 1;
  int k,l;
  double a_max;
  clock_t t0,t1;
  double dt = CLOCKS_PER_SEC;
  
  // main loop
  t0    = clock();
  a_max = getMaxOffDiag(n,A,&k,&l);
  while (true) {
    // similarity transformation
    doGivensSimTrans(n,A,V,k,l);
    a_max = getMaxOffDiag(n,A,&k,&l);
    // test whether to continue or break
    if (i>N_max || a_max<=tol) {
      t1 = clock();
      cout << "Jacobi eigenpair eigenvalue completed." << endl;
      cout << "Number of iterations = " << i << endl;
      cout << "Max. Abs. off-diag.  = " << a_max << endl;
      cout << "time spent = " << (double)(t1-t0)/dt << " s" << endl;
      break;
    };
    i++;
  };
  return V;
};

void doGivensSimTrans (const int n, double **A, double **V, int k, int l) {
  // misc
  int i;
  double s,c,tau,t;
  double a_kk = A[k][k];
  double a_ll = A[l][l];
  double a_kl = A[k][l];
  double a_ik,a_il;
  double v_ik,v_il;
  
  // compute s = sin(theta) and c = cos(theta)
  if (A[k][l] != 0.0) {
    tau = 0.5*(a_ll-a_kk)/a_kl;
    if (tau>0) {t =  1.0/( tau+sqrt(1.0+tau*tau));}
    else       {t = -1.0/(-tau+sqrt(1.0+tau*tau));};
    c = 1/sqrt(1+t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  };
  
  // compute the rotated elements
  A[k][k] = c*c*a_kk - 2.0*c*s*a_kl + s*s*a_ll;
  A[l][l] = s*s*a_kk + 2.0*c*s*a_kl + c*c*a_ll;
  A[k][l] = 0.0;
  A[l][k] = 0.0;
  
  // adjust A & V for transformation
  for (i=0; i<n; i++) {
    // elements in A other than a_kk, a_ll, a_kl, a_lk
    if (i!=k && i!=l) {
      a_ik = A[i][k];
      a_il = A[i][l];
      A[i][k] = c*a_ik - s*a_il;
      A[i][l] = c*a_il + s*a_ik;
      A[k][i] = A[i][k];
      A[l][i] = A[i][l];
    };
    // the eigenvector matrix
    v_ik = V[i][k];
    v_il = V[i][l];
    V[i][k] = c*v_ik - s*v_il;
    V[i][l] = c*v_il + s*v_ik;
  };
};


double getMaxOffDiag (const int n, double **A, int *k, int *l) {
  // misc
  int i,j;
  *k = 0;
  *l = 0;
  double max = 0.0;
  double a_abs;
  
  // test each off-diagonal element above the central diagonal
  for (i=0; i<n; i++) {
    for (j=i+1; j<n; j++) {
      a_abs = fabs(A[i][j]);
      // found a greater element
      if (a_abs > max) {
        max = a_abs;
        *k  = i;
        *l  = j;
      };
    };
  };
  return max;
};


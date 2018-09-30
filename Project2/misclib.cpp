#include "misclib.h"

/*
--------------------------------------------------------------------------------------------------------------------------
DISPLAY MATRICES & VECTORS ON THE COMMAND LINE
--------------------------------------------------------------------------------------------------------------------------
*/

// display a matrix
void display_matrix (const int m, const int n, double **M, string matrix, int precision) {
  int i,j;
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << matrix << " = " << showpos << setprecision(precision) << scientific << endl;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++) {cout << "  " << M[i][j] << "  ";};
    cout << endl;
  };
  cout.flags(f);
};

// nxm matrices
void display_matrix(const int m, const int n, double **M)                {display_matrix(m,n,M,"matrix",4);};
void display_matrix(const int m, const int n, double **M, string matrix) {display_matrix(m,n,M,matrix  ,4);};
void display_matrix(const int m, const int n, double **M, int precision) {display_matrix(m,n,M,"matrix",precision);};

// square nxn matrices
void display_matrix(const int n, double **M)                {display_matrix(n,n,M,"matrix",4);};
void display_matrix(const int n, double **M, string matrix) {display_matrix(n,n,M,matrix,  4);};
void display_matrix(const int n, double **M, int precision) {display_matrix(n,n,M,"matrix",precision);};

// display a vector
void display_vector (const int n, double *v, string vector, int precision) {
  int i,j;
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << vector << " = " << showpos << setprecision(precision) << scientific << endl;
  for (i=0; i<n; i++) {cout << "  " << v[i] << "  " << endl;};
  cout.flags(f);
};
void display_vector (const int n, double *v)                {display_vector(n,v,"vector",4);};
void display_vector (const int n, double *v, string vector) {display_vector(n,v,vector,  4);};
void display_vector (const int n, double *v, int precision) {display_vector(n,v,"vector",precision);};

/* 
--------------------------------------------------------------------------------------------------------------------------
DECLARE & DELETE MATRICES
--------------------------------------------------------------------------------------------------------------------------
*/

// initialize a zero matrix
double ** init_matrix (const int m, const int n) {
  int i,j;
  double **M;
  M = new double *[m];
  for (i=0; i<m; i++)
  {
    M[i] = new double [n];
    for (j=0; j<n; j++) {M[i][j] = 0.0;};
  };
  return M;
};
double ** init_matrix (const int n) {return init_matrix(n,n);};

// initialize a random matrix with elements between 0 and 1
double ** init_random_matrix (const int m, const int n) {
  int i,j;
  double **M;
  M = new double *[m];
  for (i=0; i<m; i++) {
    M[i] = new double [n];
    for (j=0; j<n; j++) {M[i][j] = (double) rand()/(RAND_MAX);};
  };
  return M;
};
double ** init_random_matrix (const int n) {return init_random_matrix(n,n);};

// initialize a random symmetric matrix with elements between 0 and 1
double ** init_symmetric_random_matrix (const int m, const int n) {
  int i,j;
  double **M;
  M = new double *[m];
  for (i=0; i<m; i++) {
    M[i] = new double [n];
    for (j=i; j<n; j++) {M[i][j] = (double) rand()/(RAND_MAX);};
  };
  for (i=1; i<m; i++) {for (j=0; j<i; j++) {M[i][j] = M[j][i];};};
  return M;
};
double ** init_symmetric_random_matrix (const int n) {return init_symmetric_random_matrix(n,n);};

// initialize an identity matrix
double ** init_identity_matrix (const int m, const int n) {
  int  i,l;
  double **M = init_matrix(m,n);
  if (m>=n) {l = n;}
  else      {l = m;};
  for (i=0; i<l; i++) {M[i][i] = 1;};
  return M;
};
double ** init_identity_matrix (const int n) {return init_identity_matrix(n,n);};

// initialize a tridiagonal Toeplitz matrix
double ** init_tridiag_Toeplitz_matrix(const int m, const int n, double a, double b, double c) {
  int i,j;
  double **M = init_matrix(m,n);
  M[0][0] = b;
  for (i=1; i<n; i++) {
    M[i][i-1] = a;
    M[i][i]   = b;
    M[i-1][i] = c;
  }
  return M;
};
 double ** init_tridiag_Toeplitz_matrix(const int n, double a, double b, double c) {
  return  init_tridiag_Toeplitz_matrix(n,n,a,b,c);}

// free-up space by deleting a matrix
void delete_matrix (const int n, double ** M) {
  int i;
  for(i=0; i<n; i++) {delete[] M[i];};
  delete[] M;
};

/* 
--------------------------------------------------------------------------------------------------------------------------
MATRIX & VECTOR OPERATIONS
--------------------------------------------------------------------------------------------------------------------------
*/

// scalar multiplication of each element in a matrix
void scale_matrix (const int n, double **M, double c) {
  int i,j;
  for (i=0; i<n; i++) {for (j=0; j<n; j++) {M[i][j] *= c;};};
};

// extract the diagonal from a matrix
double * get_matrix_diagonal (const int n, double **M) {
  int i;
  double *diag;
  diag = new double [n];
  for (i=0; i<n; i++) {diag[i] = M[i][i];};
  return diag;
};

// sort a vector from minimum to maximum
double * sort_vector (const int n, double *v, int *idx) {
  // misc
  int i,j;
  double c;
  double min = v[0];
  double max = v[0];
  // find max and min
  for (i=1; i<n; i++) {
    if (v[i] < min) {min=v[i];idx[0]  =i;};
    if (v[i] > max) {max=v[i];idx[n-1]=i;};
  };
  double *v_new; v_new = new double [n];
  v_new[0] = min; v_new[n-1] = max;
  // main loop - builds every element of v_new
  for (i=1; i<n-1; i++) {
    c = max;
    for (j=0; j<n; j++) {
      if (v[j]>v_new[i-1] && v[j] < c) {
        c      = v[j];
        idx[i] = j;
      };
    };
    v_new[i] = c;
  };
  return v_new;
};


#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include <string>
#include "tridiagmateq.h"

using namespace std;

double* gen_x (int n, double h) {
  int     i;
  double* x = new double[n]; 
  for (i=1; i<n+1; i++) {x[i] = i*h;}
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

int save_data (int n, string algorithm, double* x, double* exact, double* y) {
  int i;
  string   filename = to_string(n) + "_" + algorithm + ".dat";
  bool exists = file_existence("data/"+filename);
  if (exists==true) {cout << "Overwriting old file." << endl;}
  ofstream outfile("data/"+filename, ios::out | ios::trunc);
  if (!outfile.is_open()) {
    cout << "Error: Unable to open output file." << endl;
    return 1;
  }
  outfile << "n         = " << n << endl;
  outfile << "algorithm = " << algorithm << endl;
  outfile << "---------------------------------------------------------------------" << endl;
  outfile << "        x        |     y_exact     |      y_num" << endl;
  for (i=0; i<n; i++) {
    outfile << " " << left << scientific << setprecision(9) << setw(10) << x[i] << " | " << exact[i] << " | " << y[i] << endl;
  }
  outfile.close();
  return 0;
};

inline bool file_existence (const string& filename) {
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0); 
}




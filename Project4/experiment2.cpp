#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>
#include "time.h"
#include <algorithm>

using namespace std;

// period boundary conditions
inline int periodic(int i, int limit, int add) {return (i+limit+add) % (limit);};

int main (int argc, char* argv[]) {

  if (argc<5) {
    cout << "Error: Missing command line arguments" << endl;
    return 1;
  };
  
  // program parameters
  int    L = 20;        // lattice dimensions
  string outfile;       // output filename
  double T;             // temperature [kT/J]
  int    N_MC;          // number of Monte Carlo cycles
  bool   orderdisorder; // whether to start from an ordered state or a disordered state
    
  // read command line arguments
  outfile = argv[1];
  T       = atof(argv[2]);
  N_MC    = atoi(argv[3]);
  if      (string(argv[4])=="true")  {orderdisorder = true;}
  else if (string(argv[4])=="false") {orderdisorder = false;}
  else {
    cout << "Error: Invalid order/disorder argument (true/false)" << endl;
    return 1;
  };
  
  // precompute exponentials
  double * dE = new double [17];
  for (int deps=-8; deps<=8; deps+=4) {dE[deps+8] = exp(-deps/T);};
  
  // setup uniform distribution
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  // print initial message
  cout << "Running Monte Carlo experiment for " << L << " x " << L << " system with:" << endl;
  cout << "  temperature  = " << T    << endl;
  cout << "  #MC cycles   = " << N_MC << endl;
  cout << "  total #spins = " << L*L  << endl;
  if (orderdisorder==true) {cout << "  initial state: ordered"    << endl;}
  else                     {cout << "  initial state: disordered" << endl;};
  
  // init Monte Carlo integrals
  double * avg_E    = new double [N_MC+1];
  double * avg_E2   = new double [N_MC+1];
  double * abs_M    = new double [N_MC+1];
  double * avg_M    = new double [N_MC+1];
  double * avg_M2   = new double [N_MC+1];
  int *    accepted = new int    [N_MC+1];
  
  avg_E[ 0] = 0;
  avg_E2[0] = 0;
  abs_M[ 0] = 0;
  avg_M[ 0] = 0;
  avg_M2[0] = 0;

  // initialise system
  double E = 0;
  double M = 0;
  int ** state = new int * [L];
  
  if (orderdisorder==true) {
    for (int i=0; i<L; i++) {
      state[i] = new int [L];
      for (int j=0; j<L; j++) {state[i][j] = 1;};
    };
  } else {
    for (int i=0; i<L; i++) {
      state[i] = new int [L];
      for (int j=0; j<L; j++) {
        state[i][j] = ((uniform(gen)<0.5) ? (-1) : (1));
      };
    };
  };

  for (int i=0; i<L; i++) {
    for (int j=0; j<L; j++) {
      E -= (double)state[i][j] * (state[periodic(i,L,-1)][j] + state[i][periodic(j,L,-1)]);
      M += (double)state[i][j];
    };
  };
  
  // program timing - start
  clock_t start = clock();
  
  // run Metropolis sampling
  for (int cycle=0; cycle<N_MC; cycle++) {
    // sweep over every point on the lattice
    for (int x=0; x<L; x++) {
      for (int y=0; y<L; y++) {
        // propose a random spin to flip
        int i = (int) (uniform(gen)*(double)L);
        int j = (int) (uniform(gen)*(double)L);
        // compute energy difference
        int DE = 2*state[i][j] * (  state[i][periodic(j,L,-1)] + state[periodic(i,L,-1)][j]
                                  + state[i][periodic(j,L,+1)] + state[periodic(i,L,+1)][j] );
        
        accepted[cycle+1] = accepted[cycle];
        
        // accept or reject
        if (uniform(gen)<=dE[DE+8]) {
          // accept flip & update system
          state[i][j]       *= -1;
          E                 += (double) DE;
          M                 += (double) 2*state[i][j];
          accepted[cycle+1] += 1;
        };
      
      };
    };
    
    // add contribution to Monte Carlo integrals
    avg_E[ cycle+1] = avg_E[ cycle] + E;
    avg_E2[cycle+1] = avg_E2[cycle] + E*E;
    abs_M[ cycle+1] = abs_M[ cycle] + fabs(M);
    avg_M[ cycle+1] = avg_M[ cycle] + M;
    avg_M2[cycle+1] = avg_M2[cycle] + M*M;
  
  };
  
  // normalise Monte Carlo integrals
  double   norm_  = 1./((double)(L*L));
  double * norm   = new double [N_MC];
  for (int i=0; i<N_MC; i++) {
    norm[i]      = norm_/((double)(i+1));
    avg_E[ i+1] *= norm[i];
    avg_E2[i+1] *= norm[i];
    abs_M[ i+1] *= norm[i];
    avg_M[ i+1] *= norm[i];
    avg_M2[i+1] *= norm[i];
  };
  
  // program timing - end
  clock_t end     = clock();
  double  runtime = ((double)end-start)/((double)CLOCKS_PER_SEC);

  cout << "Simulation finished, writing results to file:" << endl;
  cout << "  '" << outfile << "'" << endl;
  cout << "time spent = " << runtime << " sec" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << endl;
  
  // write results to file
  ofstream File(outfile);
  File << "# Experiment 2" << endl;
  File << setprecision(16) << scientific << showpos;
  for (int i=0; i<N_MC+1; i++) {
    File << avg_E[ i]   << " ";
    File << avg_E2[i]   << " ";
    File << abs_M[ i]   << " ";
    File << avg_M[ i]   << " ";
    File << avg_M2[i]   << " ";
    File << accepted[i] << endl;
  };
  File.close();
  
  // deallocate arrays
  delete[] dE,norm,avg_E,avg_E2,abs_M,avg_M,avg_M2,accepted;
  for (int i=0; i<L; i++) {delete[] state[i];};
  delete[] state;
  
  return 0;
};



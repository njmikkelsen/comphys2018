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
  int    N_BurnIn;      // number of Monte Carlo cycles to burn in
  int    N_MC;          // number of Monte Carlo cycles to use for integration
    
  // read command line arguments
  outfile  = argv[1];
  T        = atof(argv[2]);
  N_BurnIn = atoi(argv[3]);
  N_MC     = atoi(argv[4]);
  
  // precompute exponentials
  double * dE = new double [17];
  for (int deps=-8; deps<=8; deps+=4) {dE[deps+8] = exp(-deps/T);};
  
  // setup uniform distribution
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  // print initial message
  cout << "Running Monte Carlo experiment for " << L << " x " << L << " system with:" << endl;
  cout << "  temperature     = " << T        << endl;
  cout << "  #burn-in cycles = " << N_BurnIn << endl;
  cout << "  #MC cycles      = " << N_MC     << endl;
  cout << "  total #spins    = " << L*L      << endl;
  
  // initialise Monte Carlo integrals
  double avg_E  = 0;
  double avg_E2 = 0;
  
  // initialise energy trackiing
  double * eps = new double [N_MC+1];
  eps[0]       = 0;
  
  // initialise system
  double E = 0;
  int ** state = new int * [L];  
  for (int i=0; i<L; i++) {
    state[i] = new int [L];
    for (int j=0; j<L; j++) {state[i][j] = 1;};
  };
  for (int i=0; i<L; i++) {
    for (int j=0; j<L; j++) {
      E -= (double)state[i][j] * (state[periodic(i,L,-1)][j] + state[i][periodic(j,L,-1)]);
    };
  };
  
  // run Metropolis burn-in
  cout << "Running the burn-in..." << endl;
  for (int cycle=0; cycle<N_BurnIn; cycle++) {
    // sweep over every point on the lattice
    for (int x=0; x<L; x++) {
      for (int y=0; y<L; y++) {
        // propose a random spin to flip
        int i = (int) (uniform(gen)*(double)L);
        int j = (int) (uniform(gen)*(double)L);
        // compute energy difference
        int DE = 2*state[i][j] * (  state[i][periodic(j,L,-1)] + state[periodic(i,L,-1)][j]
                                  + state[i][periodic(j,L,+1)] + state[periodic(i,L,+1)][j] );
        
        // accept or reject
        if (uniform(gen)<=dE[DE+8]) {
          // accept flip & update system
          state[i][j] *= -1;
          E           += (double) DE;
        };
      };
    };
  };
  
  cout << "Running the integration..." << endl;
  
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
        
        // accept or reject
        if (uniform(gen)<=dE[DE+8]) {
          // accept flip & update system
          state[i][j]       *= -1;
          E                 += (double) DE;
        };
      
      };
    };
    
    // add contribution to Monte Carlo integrals
    avg_E  += E;
    avg_E2 += E*E;
    
    // track energy
    eps[cycle+1] = E;
  
  };
  
  // normalise Monte Carlo integrals
  double norm1 = 1./((double)(L*L));
  double norm2 = norm1/((double)N_MC);
  for (int i=0; i<N_MC+1; i++) {eps[i] *= norm1;};
  avg_E  *= norm2;
  avg_E2 *= norm2;
  
  // program timing - end
  clock_t end     = clock();
  double  runtime = ((double)end-start)/((double)CLOCKS_PER_SEC);
  
  // compute quantities
  double eps_var = avg_E2 - avg_E*avg_E;
  double eps_std = sqrt(eps_var);
  
  // print results
  cout << "Simulation finished, writing results to file:" << endl;
  cout << "  '" << outfile << "'" << endl;
  cout << "Program output:" << endl;
  cout << "  mean energy   = " << avg_E   << endl;
  cout << "  mean energy^2 = " << avg_E2  << endl;
  cout << "  energy var    = " << eps_var << endl;
  cout << "  energy std    = " << eps_std << endl;
  cout << "  time spent    = " << runtime << " sec" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << endl;
  
  // write results to file
  ofstream File(outfile);
  File << "# Experiment 2" << endl;
  File << setprecision(16) << scientific << showpos;
  for (int i=0; i<N_MC+1; i++) {File << eps[i] << endl;};
  File.close();
  
  // deallocate arrays
  delete[] dE,avg_E,avg_E2,eps;
  for (int i=0; i<L; i++) {delete[] state[i];};
  delete[] state;
  
  return 0;
};



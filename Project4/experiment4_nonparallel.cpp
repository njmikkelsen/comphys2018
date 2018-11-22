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

  if (argc<6) {
    cout << "Error: Missing command line arguments" << endl;
    return 1;
  };
  
  // program parameters
  string outfile;     // output filename
  int      L;         // lattice dimensions
  int      N_T;       // number of temperatures
  double * T;         // temperature [kT/J]
  int      N_BurnIn;  // number of Monte Carlo cycles to burn in
  int      N_MC;      // number of Monte Carlo cycles to use for integration
  
  // read command line arguments
  outfile  = argv[1];
  L        = atoi(argv[2]);
  N_BurnIn = atoi(argv[3]);
  N_MC     = atoi(argv[4]);
  N_T      = argc-5;
  T        = new double [N_T];
  for (int i=0; i<N_T; i++) {T[i] = atof(argv[i+5]);};
  
  // setup uniform distribution
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  // open output file
  ofstream File(outfile);   
  File << "# Experiment 4 non-parallel" << endl;
  File << setprecision(16) << scientific;
  
  // print initial message
  cout << "Running " << N_T << " Monte Carlo experiments for " << L << " x " << L << " system with:" << endl;
  cout << "  #burn-in cycles = " << N_BurnIn << endl;
  cout << "  #MC cycles      = " << N_MC     << endl;
  cout << "  total #spins    = " << L*L      << endl;
  cout << "--------------------------------------------------------------------" << endl;
  
  // master loop - repeated experiments w/ different N_MC
  for (int experiment=0; experiment<N_T; experiment++) {
    
    // print initial experiment message
    cout << "experiment (" << experiment+1 << "/" << N_T << ")" << endl;
    
    // precompute exponentials
    double * dE = new double [17];
    for (int deps=-8; deps<=8; deps+=4) {dE[deps+8] = exp(-deps/T[experiment]);};
    
    // init Monte Carlo integrals
    double norm   = 1./((double)(N_MC));
    double avg_E  = 0;
    double avg_E2 = 0;
    double abs_M  = 0;
    double avg_M  = 0;
    double avg_M2 = 0;

    // initialise system
    double E = 0;
    double M = (double)L*L;
    int ** state = new int * [L];
    for (int i=0; i<L; i++) {
      state[i] = new int [L];
      for (int j=0; j<L; j++) {state[i][j]  = 1;};
    };
    for (int i=0; i<L; i++) {
      for (int j=0; j<L; j++) {
        E -= (double)state[i][j] * (state[periodic(i,L,-1)][j] + state[i][periodic(j,L,-1)]);
      };
    };

    // run Metropolis burn-in
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
            state[i][j] *= -1;
            E           += (double) DE;
            M           += (double) 2*state[i][j];
          };
        
        };
      };

      // update Monte Carlo integrals
      avg_E  += E;
      avg_E2 += E*E;
      abs_M  += fabs(M);
      avg_M  += M;
      avg_M2 += M*M;
    
    };
    
    // normalise Monte Carlo integrals
    avg_E  *= norm;
    avg_E2 *= norm;
    abs_M  *= norm;
    avg_M  *= norm;
    avg_M2 *= norm;
    
    // compute heat capacity & magnetic susceptibility
    double C_V = (avg_E2 - avg_E*avg_E)/(T[experiment]*T[experiment]);
    double chi = (avg_M2 - abs_M*abs_M)/T[experiment];
    
    // program timing - end
    clock_t end     = clock();
    double  runtime = ((double)end-start)/((double)CLOCKS_PER_SEC);
    
    // normalise per spin
    double Nspins = L*L;
    avg_E  /= Nspins;
    avg_E2 /= Nspins;
    abs_M  /= Nspins;
    avg_M  /= Nspins;
    avg_M2 /= Nspins;
    C_V    /= Nspins;
    chi    /= Nspins;
    
    // write to file
    File << avg_E   << " ";
    File << avg_E2  << " ";
    File << abs_M   << " ";
    File << avg_M   << " ";
    File << avg_M2  << " ";
    File << C_V     << " ";
    File << chi     << " ";
    File << runtime << endl;
    
    // deallocate arrays
    for (int i=0; i<L; i++) {delete[] state[i];};
    delete[] dE,state;
    
  };
  
  cout << endl;  
  
  // deallocate arrays
  delete[] T;
  
  return 0;
};


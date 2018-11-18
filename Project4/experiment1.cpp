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
  
  // program parameters
  int     L = 2;        // lattice dimensions
  double  T = 1.0;      // temperature [kT/J]
  int     N;            // number of simulations
  int *   N_MC;         // array of number of Monte Carlo cycles per simulation
  double  beta = 1./T;  // inverse temperature parameter
  string  outfile;      // path to output file (read from command line)
  
    
  // read command line arguments
  if (argc>1) {
    outfile = argv[1];
    N       = argc-2;
    N_MC    = new int [N];
    for (int i=0; i<N; i++) {
      N_MC[i] = atoi(argv[i+2]);
    };
  } else {
    cout << "Error: Missing command line arguments" << endl;
    return 1;
  };
  
  // precompute exponentials
  double * dE = new double [17];
  for (int deps=-8; deps<=8; deps+=4) {dE[deps+8] = exp(-deps*beta);};
  
  // setup uniform distribution
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  // open output file
  ofstream File(outfile);
  File << "# Experiment 1" << endl;
  File << setprecision(16) << scientific;
  int int_width = ((int)(log10(*max_element(N_MC,N_MC+N)))+1);
  
  // print initial message
  cout << "Running " << N << " Monte Carlo experiments for the " << L << " x " << L << " system with:" << endl;
  cout << "  temperature  = " << T << endl;
  cout << "  total #spins = " << L*L << endl;
  cout << "-----------------------------------------" << endl;
  cout << "Expectation values per spin:" << endl;
  
  // master loop - repeated experiments w/ different N_MC
  for (int experiment=0; experiment<N; experiment++) {
  
    // print initial experiment message
    cout << endl;
    cout << "experiment (" << experiment+1 << "/" << N << ")" << endl;
  
    // init Monte Carlo integrals
    double norm   = 1./((double)N_MC[experiment]*L*L);
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
    
    // program timing - start
    clock_t start = clock();
    
    // run Metropolis sampling
    for (int cycle=0; cycle<N_MC[experiment]; cycle++) {
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
    
    // program timing - end
    clock_t end     = clock();
    double  runtime = ((double)end-start)/((double)CLOCKS_PER_SEC);
  
    // print results
    cout << "N_MC = " << N_MC[experiment] << ":" << endl;
    cout << "  eps    = " << avg_E  << endl;
    cout << "  eps^2  = " << avg_E2 << endl;
    cout << "  |M|    = " << abs_M  << endl;
    cout << "  M      = " << avg_M  << endl;
    cout << "  M^2    = " << avg_M2 << endl;
    cout << "run-time = " << runtime << " sec" << endl;
    
    // write results to file
    File << noshowpos << setw(int_width);
    File << N_MC[experiment] << " ";
    File << showpos;
    File << avg_E   << " ";
    File << avg_E2  << " ";
    File << abs_M   << " ";
    File << avg_M   << " ";
    File << avg_M2  << " ";
    File << runtime << endl;
    
    // deallocate arrays
    for (int i=0; i<L; i++) {delete[] state[i];};
    delete[] state;
  };
  
  // close output file
  File.close();
  
  // deallocate arrays
  delete[] dE,N_MC;
  
  return 0;
};



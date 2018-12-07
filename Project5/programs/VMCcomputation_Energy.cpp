#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "vectorlib.h"
#include "vmc.h"

using namespace std;

int main (int argc, char* argv[])
{
  // verify correct number of command line arguments
  if (argc < 7) {
    cout << "Error: Missing command line arguments!" << endl;
    cout << "Expected: omega, alpha, beta, TYPE, N_MC, N_burn" << endl;
    return 1;
  };
  
  // command line variables
  double omega  = atof(argv[1]);
  double alpha  = atof(argv[2]);
  double beta   = atof(argv[3]);
  int    TYPE   = atoi(argv[4]);
  int    N_MC   = atoi(argv[5]);
  int    N_burn = atoi(argv[6]);
  
  // initialize Monte Carlo variables
  double    E_T;
  double    varE_T;
  double    acceptance;
  VMC       vmc;
  TrialWave psi(TYPE,omega);
  psi.setVarParams(alpha,beta);
    
  // run Monte Carlo simulation
  vmc.run_MonteCarlo_energy(psi,N_MC,N_burn,E_T,varE_T,acceptance);
  
  // print results
  cout << scientific << setprecision(16);
  cout << "Trial Wave Energy = " << E_T        << endl;
  cout << "Energy Variance   = " << varE_T     << endl;
  cout << "Acceptance        = " << acceptance << endl;
  
  return 0;
};

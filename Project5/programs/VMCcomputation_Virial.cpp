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
  if (argc < 8) {
    cout << "Error: Missing command line arguments!" << endl;
    cout << "Expected: interaction, omega, alpha, beta, TYPE, N_MC, N_burn" << endl;
    return 1;
  };
  
  // command line variables
  bool   interaction_;
  int    interaction = atoi(argv[1]);
  double omega       = atof(argv[2]);
  double alpha       = atof(argv[3]);
  double beta        = atof(argv[4]);
  int    TYPE        = atoi(argv[5]);
  int    N_MC        = atoi(argv[6]);
  int    N_burn      = atoi(argv[7]);
  
  // verify interaction variable
  if (interaction == 1)      {interaction_ = true;}
  else if (interaction == 0) {interaction_ = false;}
  else {
    cout << "Error: invalid interaction!"           << endl;
    cout << "Expected '0' = false or '1' = true"    << endl;
    cout << "Including the interaction by default." << endl;
    interaction_ = true;
 };
  
  // initialize Monte Carlo variables
  double    KineticEnergy;
  double    PotentialEnergy;
  double    acceptance;
  VMC       vmc;
  TrialWave psi(TYPE,omega,interaction_);
  psi.setVarParams(alpha,beta);
    
  // run Monte Carlo simulation
  vmc.run_MonteCarlo_virial(psi,N_MC,N_burn,KineticEnergy,PotentialEnergy,acceptance);
  
  // print results
  cout << scientific << setprecision(16);
  cout << "Kinetic Energy   = " << KineticEnergy   << endl;
  cout << "Potential Energy = " << PotentialEnergy << endl;
  cout << "Acceptance       = " << acceptance      << endl;
  
  return 0;
};

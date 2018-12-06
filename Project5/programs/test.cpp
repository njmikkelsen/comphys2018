#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "vectorlib.h"
#include "vmc.h"

using namespace std;

int main () {

  // system parameters
  double omega = 1.0;
  
  // trial function parameters
  double alpha = 1.0;
  
  // program constants
  double constant1 = alpha*omega;
  double constant2 = 0.5*omega*omega*(1-alpha*alpha);
  
  // Monte Carlo parameters
  int    N_MC     = 100000;
  double delta    = log(2)/constant1;
  int    accepted = 0;
  
  // initial state
  double * R0 = new double [6];
  double * R1 = new double [6];
  double r1_old  = 0;
  double r2_old  = 3*delta*delta;
  double r12_old = sqrt(r2_old);
  for (int i=0; i<3; i++) {
    R0[i]   = 0;
    R0[i+3] = delta;
  };
  
  // initiate some variables
  double r1_new;
  double r2_new;
  double r12_new;

  // initialize energy and energy^2
  double E  = 0;
  double E2 = 0;
  
  Uniform uniform1(-1.0,1.0);
  Uniform uniform2(+0.0,1.0);
  
  // Monte Carlo loop
  for (int n=0; n<N_MC; n++) {
    
    // propose new state
    for (int i=0; i<6; i++) {
      double a = uniform1();
      R1[i] = R0[i] + delta*a;;
    };
    
    // compute new state r1^2, r2^2 and r12
    double R1_30 = R1[3]-R1[0];
    double R1_41 = R1[4]-R1[1];
    double R1_52 = R1[5]-R1[2];
    
    r1_new  = R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2];
    r2_new  = R1[3]*R1[3] + R1[4]*R1[4] + R1[5]*R1[5];
    r12_new = sqrt(R1_30*R1_30 + R1_41*R1_41 + R1_52*R1_52);
    
    // compute Metropolis choice
    double alpha_M = exp(-constant1*((r1_new-r1_old) + (r2_new-r2_old)));
    
    // test proposal - update state if accepted
    double b = uniform2();
    if (b <= alpha_M) {
      accepted++;
      for (int i=0; i<6; i++) {R0[i] = R1[i];};
      r1_old  = r1_new;
      r2_old  = r2_new;
      r12_old = r12_new;
    };
    
    // compute local energy
    double E_L = 3*constant1 + constant2*(r1_old+r2_old) + 1./r12_old;
    
    // add integral contributions
    E  += E_L;
    E2 += E_L*E_L;
    
  };
  
  // normalize Monte Carlo integrals
  double E_T     = E/((double)N_MC);
  double Var_E_T = E2/((double)N_MC) - E_T*E_T;
  
  
  // print results
  cout << "Energy     = " << E_T                             << endl;
  cout << "Var[E]     = " << Var_E_T                         << endl;
  cout << "N_MC       = " << N_MC                            << endl;
  cout << "Accepted   = " << accepted                        << endl;
  cout << "Acceptance = " << ((double)accepted/(double)N_MC) << endl;
  
  
  // delete arrays
  delete[] R0,R1;
  
  return 0;
};


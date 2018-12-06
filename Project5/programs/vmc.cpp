#include "vmc.h"

/*
----------------------------------------
              VMC CLASS 
----------------------------------------
*/

// draw a random probability

double VMC::draw_prob () {return RandProb(RNG);};

// run Monte Carlo simulation

void VMC::run_MonteCarlo (TrialWave psi, int N_MC, int N_burn, double& E_T, double& varE_T, double& acceptance)
{
  // misc declerations
  double A;                           // Metropolis choice
  double b;                           // probability of acceptance
  double E_L;                         // Monte Carlo local energy
  int    a     = 0;                   // number of accepted proposals in Monte Carlo integration
  double E     = 0;                   // Monte Carlo sum: energy
  double E2    = 0;                   // Monte Carlo sum: energy^2
  double scale = 1/((double)N_MC);    // Monte Carlo normalization scale
  
  // Metropolis burn in
  for (int n=0; n<N_burn; n++) {
    // propose state and compute Metropolis choice
    psi.propose_state();
    A = psi.alphaM();
    b = draw_prob();
    // update state if proposal is accepted
    if (b <= A) {psi.update_state();};
  };
  
  // re-seed Random Number Generator (avoid running out of random numbers)
  RNG.seed((random_device())());
   
  // Monte Carlo integration
  for (int n=0; n<N_MC; n++) {
    // propose state and compute Metropolis choice
    psi.propose_state();
    A = psi.alphaM();
    b = draw_prob();
    // update state if proposal is accepted
    if (b <= A) {a++; psi.update_state();};
    // update Monte Carlo integrals
    E_L = psi.LocalEnergy();
    E  += E_L;
    E2 += E_L*E_L;
  };
  
  // normalize Monte Carlo integrals & compute acceptance
  E_T        = E*scale;
  varE_T     = E2*scale - E_T*E_T;
  acceptance = (double)a*scale;
};


/*
-------------------------------------------
            TrialWave CLASS
-------------------------------------------
*/

TrialWave::TrialWave (int type, double omega_) : omega(omega_), RNG((random_device())()) {
  // select trial wave
  if (type == 1) {
    TYPE         = type;
    alphaM_      = &TrialWave::Wave1_alphaM;
    LocalEnergy_ = &TrialWave::Wave1_LocalEnergy;
  } else if (type == 2) {
    TYPE         = type;
    alphaM_      = &TrialWave::Wave2_alphaM;
    LocalEnergy_ = &TrialWave::Wave2_LocalEnergy;
  } else {
    cout << "Error: Invalid wave TYPE. Expected 1 or 2." << endl;
    cout << "Selecting default: TYPE = 1" << endl;
    TYPE         = 1;
    alphaM_      = &TrialWave::Wave1_alphaM;
    LocalEnergy_ = &TrialWave::Wave1_LocalEnergy;
  };
  
  // initialize default variational parameters: alpha = 1, beta = 0
  setVarParams(1.0);
  
  // initialize trial wave state variables
  R1_prev.init(3);
  R1_next.init(3);
  R2_prev.init(3);
  R2_next.init(3);
  R2_prev[0] = 1.0;  // separate particles (avoids r12=0)
  
  r1_prev2 = R1_prev.norm_2();
  r2_prev2 = R2_prev.norm_2();
  r12_prev = (R2_prev-R1_prev).norm();
  
  if (TYPE == 2) {denom_prev = 1 + beta * r12_prev;};
};

// set variational parameters

void TrialWave::setVarParams (double alpha_, double beta_) {
  // adjust constants
  alpha     = alpha_;
  beta      = beta_;
  constant1 = alpha*omega;
  constant2 = 0.5*omega*omega*(1-alpha*alpha);
  // adjust random state generator
  delta = log(constant1 + 1) + 1./(constant1+0.25) - 1./1.25;
  Delta = uniform_real_distribution<double>(-delta,delta);
};
void TrialWave::setVarParams (double alpha_) {setVarParams(alpha_,0);};

// random number generation

double TrialWave::draw_scalar () {return Delta(RNG);};
Vector TrialWave::draw_vector () {
  Vector delta(3);
  for (int i=0; i<3; i++) {delta[i] = draw_scalar();};
  return delta;
};

// functions accessed by VMC method

void TrialWave::propose_state () {
  // draw random vectors
  Vector delta1 = draw_vector();
  Vector delta2 = draw_vector();
  // propose new state
  R1_next    = R1_prev + delta1;
  R2_next    = R2_prev + delta2;
  r1_next2   = R1_next.norm_2();
  r2_next2   = R2_next.norm_2();
  r12_next   = (R2_next-R1_next).norm();
  dStateNorm = r1_next2 + r2_next2 - r1_prev2 - r2_prev2;
  if (TYPE == 2) {denom_next = 1 + beta * r12_next;};
};

void TrialWave::update_state () {
  R1_prev  = R1_next;
  R2_prev  = R2_next;
  r1_prev2 = r1_next2;
  r2_prev2 = r2_next2;
  r12_prev = r12_next;
  if (TYPE == 2) {denom_prev = denom_next;};
};

double TrialWave::alphaM      () {return (this->*alphaM_)     ();};
double TrialWave::LocalEnergy () {return (this->*LocalEnergy_)();};

// wave 1 

double TrialWave::Wave1_alphaM      () {return exp( -constant1 * dStateNorm);};
double TrialWave::Wave1_LocalEnergy () {return 3*constant1 + constant2 * (r1_prev2 + r2_prev2) + 1./r12_prev;};

// wave 2

double TrialWave::Wave2_alphaM      () {return exp( -constant1 * dStateNorm + r12_next/denom_next - r12_prev/denom_prev);};
double TrialWave::Wave2_LocalEnergy () {
  double E1          = 3*constant1 + constant2 * (r1_prev2 + r2_prev2) + 1./r12_prev;
  double denom_prev2 = denom_prev*denom_prev;
  return E1 + 0.5*(constant1 * r12_prev - 2./r12_prev + 2*beta/denom_prev - 0.5/denom_prev2) / denom_prev2;
  
};


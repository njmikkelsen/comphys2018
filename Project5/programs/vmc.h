#ifndef VMC_H
#define VMC_H

#include <cstdlib>
#include <random>
#include "vectorlib.h"

using namespace std;

class VMC;
class TrialWave;

class VMC {
  private:
    // probability generator
    ranlux48 RNG;
    double   draw_prob ();
    uniform_real_distribution<double> RandProb;
    
  public:
    // constructor & destructor
    VMC  () : RNG((random_device())()), RandProb(0,1) {};
    ~VMC () {};
    
    // Monte Carlo simulations
    void run_MonteCarlo_energy     (TrialWave,int,int,double&,double&,double&);  // compute energy integrals
    void run_MonteCarlo_separation (TrialWave,int,int,double&,double&,double&);  // compute expected particle separation
    void run_MonteCarlo_virial     (TrialWave,int,int,double&,double&,double&);  // compute expected kinetic & potential energy
};

class TrialWave {
  private:
    // trial wave parameter variables
    bool   INTER;      // indicates whether the electrons should interact
    int    TYPE;       // indicates which trial wave is used, available: 1 or 2
    double omega;      // Harmonic Oscillator frequency
    double alpha;      // variational parameter
    double beta;       // variational parameter
    double constant1;  // alpha*omega
    double constant2;  // 0.5*omega*omega
    double constant3;  // 0.5*omega*omega*alpha*alpha

    // previous trial wave state variables
    Vector R1_prev;   // particle 1 position
    Vector R2_prev;   // particle 2 position
    double r1_prev2;  // norm-2 of R1_prev
    double r2_prev2;  // norm-2 of R2_prev
    double r12_prev;  // distance between particles
    
    // next trial wave state variables
    Vector R1_next;   // particle 1 position
    Vector R2_next;   // particle 2 position
    double r1_next2;  // norm-2 of R1_next    
    double r2_next2;  // norm-2 of R2_next
    double r12_next;  // distance between particles
    
    // extra trial wave state variables
    double dStateNorm;  // r1_next2 + r2_next2 - r1_prev2 - r2_prev2
    double denom_prev;  // 1 + beta*r12_prev   ('denom' as in denominator)
    double denom_next;  // 1 + beta*r12_next
    
    // random trial wave state generation
    ranlux48 RNG;                             // Random Number Generator
    uniform_real_distribution<double> Delta;  // distribution
    double delta;                             // max/min bounds of uniform distribution
    double draw_scalar ();                    // draw random number
    Vector draw_vector ();                    // draw random vector
    
    // wave-specific functions
    double Wave1_alphaM        ();  // wave 1, Metropolis choice
    double Wave2_alphaM        ();  // wave 2, Metropolis choice
    double Wave1_KineticEnergy ();  // wave 1, kinetic energy
    double Wave2_KineticEnergy ();  // wave 2, kinetic energy
    
    
    // function pointers (chosen based on TYPE and INTER)
    double (TrialWave::*alphaM_)          ();
    double (TrialWave::*LocalEnergy_)     ();
    double (TrialWave::*KineticEnergy_)   ();
    double (TrialWave::*PotentialEnergy_) ();
    
    // interaction functions
    double dont_interact ();  // returns 0
    double do_interact   ();  // returns Coulomb interaction
    
  public:
    // constructor & destructor
    TrialWave  (int,double,bool=true);
    ~TrialWave () {};
    
    // set variational parameters
    void setVarParams (double,double);    // set alpha and beta
    void setVarParams (double);           // set alpha
    
    // functions accessed by VMC class
    void   propose_state   ();
    void   update_state    ();
    double alphaM          ();
    double LocalEnergy     ();
    double KineticEnergy   ();
    double PotentialEnergy ();
    double separation      ();
};

#endif

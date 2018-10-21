#include "gravity.h"

int main () {
  
  GravitySystem SolarSystem ("The Solar System");
  
  Vector zero(3,0);
  Vector r0(3,0);
  Vector v0(3,0);
  
  // parameters
  int    N    = 10000;
  double dt   = 5e-2;
  double x0   = 1;
  double v    = 1.2e0;
//  double beta = 2.15;
  
  r0(0) = x0;
  v0(1) = v;
  
  SolarSystem.add_body("The Sun", 1.0000000e+0,zero,zero,false);
  SolarSystem.add_body("Earth",   3.0024584e-6,r0,v0);
  
  SolarSystem.relativistic_effects();
  
  SolarSystem.setup_integration(N,dt);
  SolarSystem.run("Verlet");
  
  
  SolarSystem.write_results();
  
  
  return 0;
};



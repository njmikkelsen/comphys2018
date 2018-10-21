#ifndef GRAVITY_H
#define GRAVITY_H

#include "vectorlib.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include "time.h"

class GravityBody;
class StaticBody;
class DynamicBody;
class GravitySystem;

/*
    GRAVITATIONAL BODIES
*/

class GravityBody {
  protected:
    string name_;  // name of gravitational body
    double mass_;  // mass of gravitational body
    bool   type_;  // false: static,  true: dynamic
    Vector r_;     // position vector (initial for dynamic bodies)
    Vector v_;     // velocity vector (initial for dynamic bodies)
  public:
    // constructors & destructor
    GravityBody  (bool type);
    void init    (string name, double mass, Vector r, Vector v, bool type);
    GravityBody  (string name, double mass, Vector r, Vector v, bool type);
    virtual ~GravityBody () {};
    // methods for accessing private members
    string name      ();
    double mass      ();
    string type      ();
    bool   type_bool ();
    // operator overloading: access position using () and velocity using []
    virtual Vector  operator () (int i) const = 0;  // the argument i is included
    virtual Vector& operator () (int i)       = 0;  // mainly for the dynamic bodies
    virtual Vector  operator [] (int i) const = 0;  // as it has no effect on the static
    virtual Vector& operator [] (int i)       = 0;  // bodies
    // misc
    Vector init_pos () const;
    Vector init_vel () const;
    virtual void prepare_integration (int N) {};
};

class StaticBody : public GravityBody {
  public:
    // constructors & desctructor
    StaticBody  (string name, double mass, Vector r, Vector v);
    void init   (string name, double mass, Vector r, Vector v);
    ~StaticBody ();
    // operator overloading: access position using () and velocity using []
    virtual Vector  operator () (int i) const;
    virtual Vector& operator () (int i);
    virtual Vector  operator [] (int i) const;
    virtual Vector& operator [] (int i);
};

class DynamicBody : public GravityBody {
  private:
    TimeVector R_;  // time-dependent position vector
    TimeVector V_;  // time-dependent velocity vector
    int        N_;  // number of time points (length of R_ and V_)
  public:
    //constructors & destructor
    DynamicBody  (string name, double mass, Vector r, Vector v);
    void init    (string name, double mass, Vector r, Vector v);
    ~DynamicBody ();
    // operator overloading: access position using () and velocity using []
    virtual Vector  operator () (int i) const;
    virtual Vector& operator () (int i);
    virtual Vector  operator [] (int i) const;
    virtual Vector& operator [] (int i);
    // misc
    void prepare_integration (int N);
};

/*
    GRAVITATIONAL SYSTEMS
*/

class GravitySystem {
  private:
    // variables
    double G;   // the gravitational constant in SI units
    double c;   // the speed of light         in SI units
    double au;  // the astronomical unit      in SI units
    double yr;  // the Julian year            in SI units

    string  sysname_;  // name of gravitational system
    double  G_;        // effective gravitational constant
    double  c_;        // effective light speed
    bool    run_;      // denotes whether the simulation has been run
    bool    ready_;    // denotes whether the initial conditions are ready
    int     N_tot;     // total number of gravitational bodies
    int     N_dyn;     // number of dynamic bodies in the system
    int     N_sta;     // number of static  bodies in the system
    int     N_int;     // number of integration steps (length of R/V + 1)
    double  dt_;       // step length of integration
    clock_t t_0;       // time at start of simulation
    clock_t t_1;       // time at end of simulation
    double  beta_;     // constant in gravitational acceleration proportionality
    
    Vector* g_Verlet;  // past gravity vector used in the velocity Verlet method
    
    vector<int>          idx_dyn;  // vector of indices for the dynamic bodies
    vector<int>          idx_sta;  // vector of indices for the static bodies
    vector<GravityBody*> Bodies;   // vector of pointers to gravitational bodies in the system
    // methods
    int find_body_idx (string name); // locate index of body "name" in Bodies, returns -1 if "name" not in Bodies

    Vector* positions  (int i);  // returns the positions  of every body at time index i
    Vector* velocities (int i);  // returns the velocities of every body at time index i

    Vector* (GravitySystem::*gravity) (int i, Vector* POS);  // pointer to gravitational acceleration
    Vector* beta_gravity              (int i, Vector* POS);  // returns Newtonian-like gravity proportional to 1/r^beta, default beta=2
    Vector* relativistic_gravity      (int i, Vector* POS);  // returns Newtonian gravity with a relativistic correction
    
    void (GravitySystem::*step)(int i);  // pointer to integration step method: Euler or Verlet
    void Euler                 (int i);  // Euler's method
    void Verlet                (int i);  // the velocity Verlet method
  public:
    // constructors and destructor
    GravitySystem  (string sysname, bool SI_units=false);
    ~GravitySystem ();
    // choose units
    void SI ();
    void AU ();
    // methods for accessing private members
    string sysname       ();
    double gravity_const ();
    double light_speed   ();
    int    N             ();
    int    N_dynamic     ();
    int    N_static      ();
    void   print_info (int precision=4);
    // methods for adjusting the bodies in the system
    void add_body      (string name, double mass, Vector r0, Vector v0, bool dynamic=true);
    void remove_body   (string name);
    void switch_status (string name);
    // methods connected to the integration
    void adjust_gravity       (double beta);
    void relativistic_effects ();
    void setup_integration    (int N, double dt);
    void run                  (string integration="Verlet");
    // misc methods
    void write_results (int precision=12);
};

#endif

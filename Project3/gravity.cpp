#include "gravity.h"

/*
-------------------------------------
          GravityBody Class
--------------------------------------
*/

GravityBody::GravityBody (bool type) {
  name_ = "";
  mass_ = 0;
  type_ = type;
  r_.init(3,0);
  v_.init(3,0);
};

void GravityBody::init (string name, double mass, Vector r, Vector v, bool type) {
  name_ = name;
  mass_ = mass;
  r_.init(3,0);
  v_.init(3,0);
  r_    = r;
  v_    = v;
  type_ = type;
};

GravityBody::GravityBody (string name, double mass, Vector r, Vector v, bool type) {
  name_ = name;
  mass_ = mass;
  r_.init(3,0);
  v_.init(3,0);
  r_    = r;
  v_    = v;
  type_ = type;
};

string GravityBody::name      () {return name_;};
double GravityBody::mass      () {return mass_;};
string GravityBody::type      () {return (type_) ? ("dynamic") : ("static");};
bool   GravityBody::type_bool () {return type_;};

Vector GravityBody::init_pos () const {return r_;};
Vector GravityBody::init_vel () const {return v_;};


/*
-------------------------------------
          StaticBody Class
--------------------------------------
*/

StaticBody::StaticBody (string name, double mass, Vector r, Vector v) : GravityBody(name,mass,r,v,false) {};

void StaticBody::init (string name, double mass, Vector r, Vector v) {GravityBody::init(name,mass,r,v,false);};

StaticBody::~StaticBody () {
  name_ = "";
  mass_ = 0;
};

Vector  StaticBody::operator () (int i) const {return r_;};
Vector& StaticBody::operator () (int i)       {return r_;};
Vector  StaticBody::operator [] (int i) const {return v_;};
Vector& StaticBody::operator [] (int i)       {return v_;};



/*
-------------------------------------
         DynamicBody Class
--------------------------------------
*/

DynamicBody::DynamicBody (string name, double mass, Vector r, Vector v) : GravityBody(name,mass,r,v,true) {};

void DynamicBody::init (string name, double mass, Vector r, Vector v) {GravityBody::init(name,mass,r,v,true);};

DynamicBody::~DynamicBody () {
  name_ = "";
  mass_ = 0;
  N_    = 0;
};

Vector  DynamicBody::operator () (int i) const {return R_(i);};
Vector& DynamicBody::operator () (int i)       {return R_(i);};
Vector  DynamicBody::operator [] (int i) const {return V_(i);};
Vector& DynamicBody::operator [] (int i)       {return V_(i);};

void DynamicBody::prepare_integration (int N) {
  N_ = N;
  R_.init(3,N);
  V_.init(3,N);
  R_(0) = r_;
  V_(0) = v_;
};


/*
-------------------------------------
        GravitySystem Class
--------------------------------------
*/

GravitySystem::GravitySystem (string sysname, bool SI_units) {
  G  = 6.67408e-11;      // m^3 kg^-1 s^-2
  c  = 299792458.;       // m s^-1
  au = 149597870700.;    // m
  yr = 365.25*24*60*60;  // s
  if (SI_units) {SI();}
  else          {AU();};
  sysname_ = sysname;
  run_     = false;
  ready_   = false;
  N_tot    = 0;
  N_dyn    = 0;
  N_sta    = 0;
  N_int    = 0;
  dt_      = 0;
  beta_    = 2;
  gravity  = &GravitySystem::beta_gravity;
};

void GravitySystem::SI () {
  G_ = G;
  c_ = c;
};

void GravitySystem::AU () {
  G_  = 39.47841760435743;  // 4pi^2
  c_  = c*yr/au;
};

GravitySystem::~GravitySystem () {
  G        = 0;
  c        = 0;
  au       = 0;
  yr       = 0;
  G_       = 0;
  c_       = 0;
  sysname_ = "";
  run_     = false;
  ready_   = false;
  N_tot    = 0;
  N_dyn    = 0;
  N_sta    = 0;
  N_int    = 0;
  dt_      = 0;
  for (int i=0; i<N_tot; i++) {delete Bodies[i];};
  Bodies.clear();
};

string GravitySystem::sysname       () {return sysname_;};
double GravitySystem::gravity_const () {return G_;};
double GravitySystem::light_speed   () {return c_;};
int    GravitySystem::N             () {return N_tot;};
int    GravitySystem::N_dynamic     () {return N_dyn;};
int    GravitySystem::N_static      () {return N_sta;};

void GravitySystem::print_info (int precision) {
  ios_base::fmtflags f(cout.flags());
  cout.flags(f);
  cout << showpos << setprecision(precision)        << scientific;
  cout << "Gravity system: '"           << sysname_ << "':" << endl;
  cout << "Gravitational constant   = " << G_       <<         endl << noshowpos;
  cout << "total number of bodies   = " << N_tot    <<         endl;
  cout << "number of dynamic bodies = " << N_dyn    <<         endl;
  cout << "number of static  bodies = " << N_sta    <<         endl;
  cout.flags(f);
};

int GravitySystem::find_body_idx (string name) {
  int k=-1;
  for (int i=0; i<N_tot; i++) {
    if (Bodies[i]->name()==name) {k=i;};
  };
  return k;
};

void GravitySystem::add_body (string name, double mass, Vector r, Vector v, bool dynamic) {
  // check if body already added
  if (find_body_idx(name)>=0) {
    cout << "Body '" << name << "' has already been added to the system!" << endl;
  } else {
    // verify mass is positive
    if (mass<=0) {
      cout << "Body '" << name << "' cannot have negative mass!" << endl;
    } else {
      // add dynamic/static body to Bodies
      if (dynamic) {
        DynamicBody* body = new DynamicBody(name,mass,r,v);
        Bodies.push_back(body);
        idx_dyn.push_back(N_tot);
        N_dyn += 1;
      } else {
        Bodies.push_back(new StaticBody(name,mass,r,v));
        idx_sta.push_back(N_tot);
        N_sta += 1;
      };
      N_tot += 1;
    };
  };
};

void GravitySystem::remove_body (string name) {
  // check if name in Bodies
  int k = find_body_idx(name);
  if (k>=0) {
    // remove body
    if (Bodies[k]->type_bool()==true) {
      idx_dyn.erase(remove(idx_dyn.begin(),idx_dyn.end(),k),idx_dyn.end());
      N_dyn -= 1;
    } else {
      idx_sta.erase(remove(idx_sta.begin(),idx_sta.end(),k),idx_sta.end());
      N_sta -= 1;
    };
    Bodies.erase(Bodies.begin()+k);
    N_tot -= 1;
  } else {
    cout << "Body '" << name << "' does not exist in this system!" << endl;
  };
};

void GravitySystem::switch_status (string name) {
  // check if name in Bodies
  int k = find_body_idx(name);
  if (k>=0) {
    double m  = Bodies[k]->mass();
    bool   t  = Bodies[k]->type_bool();
    Vector rr = Bodies[k]->init_pos();
    Vector vv = Bodies[k]->init_vel();
    remove_body(name);
    add_body(name,m,rr,vv,!t);
  } else {
    cout << "Body '" << name << "' does not exist in this system!" << endl;
  };
};

void GravitySystem::adjust_gravity (double beta) {
  gravity = &GravitySystem::beta_gravity;
  beta_   = beta;
};

void GravitySystem::relativistic_effects () {
  gravity = &GravitySystem::relativistic_gravity;
  beta_   = 2;
};

void GravitySystem::setup_integration (int N, double dt) {
  ready_ = true;
  N_int  = N;
  dt_    = dt;
  for (int i=0; i<N_dyn; i++) {
    int k = idx_dyn[i];
    Bodies[k]->prepare_integration(N+1);
  };
};

void GravitySystem::run (string integration) {
  if (ready_) {
    // choose integration method
    if (integration=="Euler") {
      step = &GravitySystem::Euler;
    } else {
      if (integration!="Verlet") {cout << "Invalid integration scheme! Using 'Verlet'." << endl;};
      step        = &GravitySystem::Verlet;
      Vector* POS = positions(0);
      g_Verlet    = (this->*gravity)(0,POS);
    };
    // integration
    t_0 = clock();
    for (int i=0; i<N_int; i++) {(this->*step)(i);};
    t_1 = clock();
    run_ = true;
  } else {
    cout << "Error: System is not ready! Missing integration parameters." << endl;
  };
};

Vector* GravitySystem::positions (int i) {
  Vector* POS = new Vector [N_tot];
  for (int j=0; j<N_tot; j++) {
    POS[j].init(3);
    POS[j] = (*Bodies[j])(i);
  };
  return POS;
};

Vector* GravitySystem::velocities (int i) {
  Vector* VEL = new Vector [N_tot];
  for (int j=0; j<N_tot; j++) {
    VEL[j].init(3);
    VEL[j] = (*Bodies[j])[i];
  };
  return VEL;
};

Vector* GravitySystem::beta_gravity (int i, Vector* POS) {
  Vector* g = new Vector [N_dyn];
  Vector dr(3,0);
  double dr_ = 0;
  for (int j=0; j<N_dyn; j++) {   // cycle through each dynamic body
    int l = idx_dyn[j];
    g[j].init(3);
    for (int k=0; k<N_tot; k++) {   // cycle through all bodies
      if (k!=l) {
        dr    = POS[k]-POS[l];
        dr_   = dr.norm();
        g[j] += (*Bodies[k]).mass()*dr/pow(dr_,beta_+1);
      };
    };
    g[j] *= G_;
  };
  return g;
};

Vector* GravitySystem::relativistic_gravity (int i, Vector* POS) {
  Vector* g   = beta_gravity(i,POS);
  Vector* VEL = velocities(i);
  for (int j=0; j<N_dyn; j++) {
    int k  = idx_dyn[j];
    g[j]  *= 1 + 3*(POS[k]%VEL[k]).norm_2()/(POS[k].norm_2()*c_*c_);
  };
  return g;
};

void GravitySystem::Euler (int i) {
  Vector* POS = positions(i);
  Vector* g   = (this->*gravity)(i,POS);
  for (int j=0; j<N_dyn; j++) {
    int k             = idx_dyn[j];
    (*Bodies[k])(i+1) = (*Bodies[k])(i) + dt_*(*Bodies[k])[i];  // update position
    (*Bodies[k])[i+1] = (*Bodies[k])[i] + dt_*g[j];             // update velocity
  };
};

void GravitySystem::Verlet (int i) {
  // update position
  for (int j=0; j<N_dyn; j++) {
    int k             = idx_dyn[j];
    (*Bodies[k])(i+1) = (*Bodies[k])(i) + dt_*(*Bodies[k])[i] + 0.5*dt_*dt_*g_Verlet[j];
  };
  // update velocity
  Vector* POS = positions(i+1);
  Vector* g   = (this->*gravity)(i,POS);
  for (int j=0; j<N_dyn; j++) {
    int k             = idx_dyn[j];
    (*Bodies[k])[i+1] = (*Bodies[k])[i] + 0.5*dt_*(g_Verlet[j]+g[j]);
    g_Verlet[j]       = g[j];
  };
};

void GravitySystem::write_results (int precision) {
  // create/overwrite folder
  int success = system(("python3 prepare_folder.py '" + sysname_+"'").c_str());
  if (success < 0) {
    cout << "Error! Filed to run python script 'prepare_folder.py'" << endl;
  } else {
    string dirpath = "./results/" + sysname_ + "/";
    // write general results
    ofstream   File_g(dirpath + "general.dat");
    streambuf* coutbuf = cout.rdbuf();  // save cout buf
    cout.rdbuf(File_g.rdbuf());         // redirect cout stream to File
    print_info(precision);              // run print_info()
    cout.rdbuf(coutbuf);                // reset cout stream
    File_g << "no. of time units spent  = " << t_1-t_0        << endl;
    File_g << "length of time unit      = " << CLOCKS_PER_SEC << endl;
    File_g << setprecision(precision) << scientific;
    File_g << "no. of integration steps = " << N_int          << endl;
    File_g << "time step length         = " << dt_            << endl;
    File_g.close();
    // write results for each dynamic body
    for (int i=0; i<N_dyn; i++) {
      int k = idx_dyn[i];
      ofstream File(dirpath + Bodies[k]->name() + ".dat");
      File << showpos << setprecision(precision) << scientific;
      for (int j=0; j<N_int+1; j++) {
        File << (*Bodies[k])(j)(0) << "  " << (*Bodies[k])(j)(1) << "  "  << (*Bodies[k])(j)(2) << "  ";
        File << (*Bodies[k])[j](0) << "  " << (*Bodies[k])[j](1) << "  "  << (*Bodies[k])[j](2) << endl;
      };
      File.close();
    };
    // write results for static bodies
    ofstream File_s(dirpath + "static_bodies.dat");
    File_s << showpos << setprecision(precision) << scientific;
    for (int i=0; i<N_sta; i++) {
      int k = idx_sta[i];
      File_s << Bodies[k]->name() << ":" << endl;
      File_s << (*Bodies[k])(0)(0) << "  " << (*Bodies[k])(0)(1) << "  "  << (*Bodies[k])(0)(2) << "  ";
      File_s << (*Bodies[k])[0](0) << "  " << (*Bodies[k])[0](1) << "  "  << (*Bodies[k])[0](2) << endl;
    };
    File_s.close();
  };
};


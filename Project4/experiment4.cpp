#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <random>
#include <mpi.h>

using namespace std;

// period boundary conditions
inline int periodic(int i, int limit, int add) {return (i+limit+add) % (limit);};

// load comand line arguments to program variables
int load_command_line_arguments (int, char**, string&, ofstream&, int&, int&, int&, int&, double&, double&, double&, int&);

// compute the ground state energy and magnetisation
void compute_ground_state (const int, double&, double&);

// initialises and returns an Ising ground state
int ** init_ground_state (const int L);

// performs a single Metropolis sampling (single step in Markov chain)
void Metropolis_Sampling (const int Nspins, const int L, int ** state, double& E, double& M, double * dE,
                          uniform_real_distribution<double> uniform, mt19937_64 gen); 

// main program
int main (int argc, char* argv[])
{
  // program parameters
  string   outfile;  // name of results file
  ofstream File;     // output file
  double   t0;       // time at MPI computation init
  double   t1;       // time at MPI computation end
  double   Dt;       // time taken by MPI computation
  
  // Ising parameters
  int    L;       // square lattice side length
  int    Nspins;  // number of spins in lattice (=L*L)
  double E0;      // ground state energy
  double M0;      // ground state magnetisation
  double T0;      // minimum temperature
  double T1;      // maximum temperature
  double DT;      // size of temperature range
  int    N_T;     // number of temperatures considered
  
  // Metropolis parameters
  int N_burnin;  // number of burn-in cycles
  int N_MC;      // nubmer of Monte Carlo cycles
  
  // MPI parameters
  int N_procs;  // number of processes
  int rank;     // local process rank

  // init MPI computation
  MPI_Init      (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &N_procs);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  t0 = MPI_Wtime();
  
  // setup random number generator
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> uniform(0.0,1.0);
  
  // load command line arguments
  int error = load_command_line_arguments(argc,argv,outfile,File,L,Nspins,N_burnin,N_MC,T0,T1,DT,N_T);
  if (error==1) {return 1;};
  
  // setup the ground state energy & magnetistion
  compute_ground_state(L,E0,M0);
  
  // local MPI variables
  int    N_T_local = N_T/N_procs;
  double T0_local  = T0       +  rank*N_T_local*DT;
  double T1_local  = T0_local +       N_T_local*DT;
  double T_local   = T0_local;
  
  // loop over local temperature range
  for (int n=0; n<N_T_local; n++) {
  
    // precompute local exponentials
    double * dE = new double [17];
    for (int deps=-8; deps<=8; deps+=4) {dE[deps+8] = exp(-deps/T_local);};
    
    // init Monte Carlo integrals
    double avg_E  = 0;
    double avg_E2 = 0;
    double abs_M  = 0;
    double avg_M  = 0;
    double avg_M2 = 0;
    
    // initialise state
    double E = E0;
    double M = M0;
    int ** state = init_ground_state(L);
    
    // burn in
    for (int cycle=0; cycle<N_burnin; cycle++) {
      Metropolis_Sampling(Nspins,L,state,E,M,dE,uniform,gen);
    };
    
    // Monte Carlo integration
    for (int cycle=0; cycle<N_MC; cycle++) {
      // Markov chain step
      Metropolis_Sampling(Nspins,L,state,E,M,dE,uniform,gen);
      
      // update Monte Carlo integrals
      avg_E  += E;
      avg_E2 += E*E;
      abs_M  += fabs(M);
      avg_M  += M;
      avg_M2 += M*M;
    };
    
    // normalise Monte Carlo integrals
    avg_E  /= ((double)(N_MC));
    avg_E2 /= ((double)(N_MC));
    abs_M  /= ((double)(N_MC));
    avg_M  /= ((double)(N_MC));
    avg_M2 /= ((double)(N_MC));
    
    // compute heat capacity & magnetic susceptibility
    double C_V = (avg_E2 - avg_E*avg_E)/(T_local*T_local);
    double chi = (avg_M2 - abs_M*abs_M)/T_local;
    
    // normalise per spin
    avg_E  /= (double)(Nspins);
    avg_E2 /= (double)(Nspins);
    abs_M  /= (double)(Nspins);
    avg_M  /= (double)(Nspins);
    avg_M2 /= (double)(Nspins);
    C_V    /= (double)(Nspins);
    chi    /= (double)(Nspins);
    
    // write to file
    File << T_local << " ";
    File << avg_E   << " ";
    File << avg_E2  << " ";
    File << abs_M   << " ";
    File << avg_M   << " ";
    File << avg_M2  << " ";
    File << C_V     << " ";
    File << chi     << endl;
    
    // deallocate arrays
    for (int i=0; i<L; i++) {delete[] state[i];};
    delete[] dE,state;
    
    cout << "This is rank " << rank << ". I have completed T = " << T_local << endl;
    
    // update temperature
    T_local += DT;
    
  };
  
  // wait for all ranks to be finished
  MPI_Barrier (MPI_COMM_WORLD);
  
  // finish program
  t1 = MPI_Wtime();
  Dt = t1-t0;
  
  if (rank==0) {
    cout << "Computation finished." << endl;
    cout << "This program spent " << Dt << " seconds with " << N_procs << " number of processes." << endl;
  }
  
  // End MPI computation
  MPI_Finalize ();  
  
  return 0;
};



/*
--------------------------
    PROGRAM FUNCTIONS
--------------------------
*/


int load_command_line_arguments (int argc, char** argv, string& outfile, ofstream& File,
                                  int& L, int& Nspins, int& N_burnin, int& N_MC, double& T0, double& T1, double& DT, int& N_T)
{

  // verify necessary number of arguments
  if (argc<8)
  {
    cout << "Error: Missing command line arguments" << endl;
    cout << "Expected arguments:" << endl;
    cout << "prog.x output_filename lattice_size N_BurnIn N_MC T0 T1 N_T" << endl;
    return 1;
  };  
  
  // load command line arguments
  outfile  = argv[1];
  L        = atoi(argv[2]);
  N_burnin = atoi(argv[3]);
  N_MC     = atoi(argv[4]);
  T0       = atof(argv[5]);
  T1       = atof(argv[6]);
  N_T      = atoi(argv[7]);
  
  // compute resulting variables
  Nspins = L*L;
  DT     = (T1-T0)/((double)(N_T-1));
  
  // open output file
  File.open(outfile, ofstream::out | ofstream::app);
  File << setprecision(16) << scientific;
  
  return 0;
};


void compute_ground_state (const int L, double& E0, double& M0)
{
  // init
  M0 = (double)(L*L);
  E0 = 0;
  int ** state = init_ground_state(L);
  // compute energy
  for (int i=0; i<L; i++) {
    for (int j=0; j<L; j++) {
      E0 -= (double)state[i][j] * (state[periodic(i,L,-1)][j] + state[i][periodic(j,L,-1)]);
    };
  };
  // deallocate
  for (int i=0; i<L; i++) {delete[] state[i];};
  delete[] state;
};


int ** init_ground_state (const int L)
{
  int ** state = new int * [L];
  for (int i=0; i<L; i++) {
    state[i] = new int [L];
    for (int j=0; j<L; j++) {state[i][j]  = 1;};
  };
  return state;
};


void Metropolis_Sampling (const int Nspins, const int L, int ** state, double& E, double& M, double * dE,
                          uniform_real_distribution<double> uniform, mt19937_64 gen)
{
  
  // draw L*L random spins
  for (int n=0; n<Nspins; n++) {
  
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



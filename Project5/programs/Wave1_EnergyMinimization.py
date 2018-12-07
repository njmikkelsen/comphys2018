import sys,os
from subprocess import PIPE, run
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from progress_bar import progress_bar

# adjust saved figure dpi
matplotlib.rcParams.update({'savefig.dpi':300})

# make C++ executable
def make(program):
  os.system("make PROG=VMCcomputation_{:s}.o".format(program))

# run program and extract command line output
def run_program(omega,alpha,beta,N_MC,N_burn):
  result = run("./run.x {:f} {:f} {:f} 1 {:d} {:d}".format(omega,alpha,beta,N_MC,N_burn),
               stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
  return [float(line.split("=")[-1]) for line in list(filter(None,result.stdout.split("\n")))]

"""
Program for performing a Variational Monte Carlo analysis for two Coulomb-interacting
electrons in a harmonic oscillator. The VMC is done using the C++ programs:
"VMCcomputation_Energy.cpp" and "VMCcomputation_Separation.cpp".

Variational wave function:

  psi1 = exp( - 0.5 * alpha * omega * ( (r_1)^2 + (r_2)^2 ) )
  
where (r_i)^2 is the norm-2 of the position vector of particle i.

The C++ programs' command line arguments:  ./run.x omega alpha beta TYPE N_MC N_burn
where:
  omega  | harmonic oscillator frequency
  alpha  | variational parameter
  beta   | variational parameter (used with psi2, see Wave2_EnergyMinimization.py)
  TYPE   | indicates whether psi1 or psi2 is used (=1)
  N_MC   | number of Monte Carlo cycles
  N_burn | number of Metropolis burn-in cycles
"""

# quick-selection
run_section1 = False   # Analyse Monte Carlo stability
run_section2 = False   # Analyse the variational parameter space w/ respect to energy
run_section3 = True    # Analyse expected separation between the electrons

###
### Section 1: Analyse Monte Carlo stability
###

if run_section1:  
  # data set parameters
  tag          = "1"     # extra save tag used in np.save/np.load of results
  load_results = True    # if True, skip run computations
  
  # load results or run computations
  if load_results:
    N_MC,TrialEnergy,EnergyVariance,Acceptance = np.load("../results/wave1/stability_{:s}.npy".format(tag))
  else:
    # analysis parameters
    repetitions  = 20      # number of times the VMC is run with each N_MC value
    MinExponent  = 2       # 10^MinExponent is the smallest N_MC value
    MaxExponent  = 7       # 10^MaxExponent is the largest N_MC value
    numN_MC      = 20      # number of different N_MC values
    
    # misc parameters
    omega  = 1.0  # harmonic oscillator frequency
    alpha  = 1.0  # variational parameter
    beta   = 0.0  # variational parameter
    N_burn = 0    # Metropolis burn-in
    
    # logarithmically spaced array of N_MC values
    N_MC = np.logspace(2,MaxExponent,numN_MC).astype(int)
    
    # result arrays
    TrialEnergy    = np.zeros((numN_MC,repetitions))
    EnergyVariance = np.zeros((numN_MC,repetitions))
    Acceptance     = np.zeros((numN_MC,repetitions))
    
    # run computations
    make("Energy")
    sanity = progress_bar(numN_MC*repetitions)
    for j in range(repetitions):
      for i,n_MC in enumerate(N_MC):
        TrialEnergy[i,j],EnergyVariance[i,j],Acceptance[i,j] = run_program(1,1,0,n_MC,0)  # omega = alpha = 1, beta = 0, N_burn = 0
        sanity.update()
    N_MC = np.tile(N_MC,(repetitions,1)).T
    np.save("../results/wave1/stability_{:s}.npy".format(tag),[N_MC,TrialEnergy,EnergyVariance,Acceptance])
  
  # scatter-plot TrialEnergy, EnergyVariance and acceptance
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,8))
  
  fig.suptitle("Variational Monte Carlo Stability")
  
  axes[0].scatter(N_MC,TrialEnergy,   s=2.0,c='r',marker='o')
  axes[1].scatter(N_MC,EnergyVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(N_MC,Acceptance,    s=2.0,c='r',marker='o')
  
  axes[0].set_title("trial energy")
  axes[1].set_title("energy variance")
  axes[2].set_title("acceptance")
  
  axes[0].set_ylabel("trial energy")
  axes[1].set_ylabel("energy variance")
  axes[2].set_ylabel("acceptance")
  axes[2].set_xlabel("Monte Carlo cycles")
  
  axes[2].set_xlim([N_MC.min()*0.9,N_MC.max()*1.1])
  
  axes[0].set_xscale("log")
  axes[1].set_xscale("log")
  axes[2].set_xscale("log")
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.9)
  
  fig.savefig("../results/wave1/stability_{:s}.png".format(tag))
  plt.show()


###
### Section 2: Analyse the variational parameter space w/ respect to energy
###

if run_section2:
  # data set parameters
  tag          = "1"     # extra save tag used in np.save/np.load of results
  load_results = True    # if True, skip run computations
  
  # load results or run computations
  if load_results:
    alpha,TrialEnergy,EnergyVariance,Acceptance = np.load("../results/wave1/energy_minimization_{:s}.npy".format(tag))
  else:
    # parameters based on Section 1:
    N_burn       = 10000   # Metropolis burn-in
    N_MC         = 10000   # number of Monte Carlo cycles
    
    # misc parameters
    omega = 1.0  # harmonic oscillator frequency
    beta  = 0.0  # variational parameter
    
    # analysis parameters
    alpha_min   = 0.5   # minimum alpha
    alpha_max   = 1.5   # maximum alpha
    N_alpha     = 100   # number of alpha parameters
    repetitions = 5     # number of times the VMC is run with each alpha value
    
    # logarithmically spaced array of alpha values
    alpha = np.linspace(alpha_min,alpha_max,N_alpha)
    
    # result arrays
    TrialEnergy    = np.zeros((N_alpha,repetitions))
    EnergyVariance = np.zeros((N_alpha,repetitions))
    Acceptance     = np.zeros((N_alpha,repetitions))
    
    # run computations
    make("Energy")
    sanity = progress_bar(N_alpha*repetitions)
    for j in range(repetitions):
      for i,a in enumerate(alpha):
        TrialEnergy[i,j],EnergyVariance[i,j],Acceptance[i,j] = run_program(omega,a,beta,N_MC,N_burn)
        sanity.update()
    alpha = np.tile(alpha,(repetitions,1)).T
    np.save("../results/wave1/energy_minimization_{:s}.npy".format(tag),[alpha,TrialEnergy,EnergyVariance,Acceptance])
  
  # scatter-plot TrialEnergy, EnergyVariance and acceptance
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,8))
  
  fig.suptitle("Variational Monte Carlo Energy Minimization")
  
  axes[0].scatter(alpha,TrialEnergy,   s=2.0,c='r',marker='o')
  axes[1].scatter(alpha,EnergyVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(alpha,Acceptance,    s=2.0,c='r',marker='o')
  
  axes[0].set_title("trial energy")
  axes[1].set_title("energy variance")
  axes[2].set_title("acceptance")
  
  axes[0].set_ylabel("trial energy")
  axes[1].set_ylabel("energy variance")
  axes[2].set_ylabel("acceptance")
  axes[2].set_xlabel(r"$\alpha$")
  
  axes[0].set_yscale("log")
  axes[1].set_yscale("log")
  
  axes[2].set_xlim([alpha.min()*0.9,alpha.max()*1.1])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.9)
  
  fig.savefig("../results/wave1/energy_minimization_{:s}.png".format(tag))
  plt.show()

###
### Section 3: Analyse expected separation between the electrons
###

if run_section3:
  # data set parameters
  tag          = "1"     # extra save tag used in np.save/np.load of results
  load_results = True    # if True, skip run computations
  
  # load results or run computations
  if load_results:
    omega,Separation,SeparationVariance,Acceptance = np.load("../results/wave1/expected_separation_{:s}.npy".format(tag))
  else:
    # parameters based on Sections 1 and 2:
    N_burn = 10000   # Metropolis burn-in
    N_MC   = 10000   # number of Monte Carlo cycles
    alpha  = 0.85    # variational parameter
    
    # misc parameters
    beta  = 0.0  # variational parameter
    
    # analysis parameters
    omega_min   = 0.05   # minimum alpha
    omega_max   = 1.00   # maximum alpha
    N_omega     = 100    # number of alpha parameters
    repetitions = 5      # number of times the VMC is run with each omega value
    
    # logarithmically spaced array of alpha values
    omega = np.linspace(omega_min,omega_max,N_omega)
    
    # result arrays
    Separation         = np.zeros((N_omega,repetitions))
    SeparationVariance = np.zeros((N_omega,repetitions))
    Acceptance         = np.zeros((N_omega,repetitions))
    
    # run computations
    make("Separation")
    sanity = progress_bar(N_omega*repetitions)
    for j in range(repetitions):
      for i,om in enumerate(omega):
        Separation[i,j],SeparationVariance[i,j],Acceptance[i,j] = run_program(om,alpha,beta,N_MC,N_burn)
        sanity.update()
    omega = np.tile(omega,(repetitions,1)).T
    np.save("../results/wave1/expected_separation_{:s}.npy".format(tag),[omega,Separation,SeparationVariance,Acceptance])
  
  # scatter-plot TrialEnergy, EnergyVariance and acceptance
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,8))
  
  fig.suptitle("Variational Monte Carlo Expected Separation Between Electrons")
  
  axes[0].scatter(omega,Separation,        s=2.0,c='r',marker='o')
  axes[1].scatter(omega,SeparationVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(omega,Acceptance,        s=2.0,c='r',marker='o')
  
  axes[0].set_title("expected separation")
  axes[1].set_title("separation variance")
  axes[2].set_title("acceptance")
  
  axes[0].set_ylabel("expected separation")
  axes[1].set_ylabel("separation variance")
  axes[2].set_ylabel("acceptance")
  axes[2].set_xlabel(r"$\omega$")
  
  axes[2].set_xlim([omega.min()*0.9,omega.max()*1.1])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.9)
  
  fig.savefig("../results/wave1/expected_separation_{:s}.png".format(tag))
  plt.show()

# clean up .o files
os.system("make clean")

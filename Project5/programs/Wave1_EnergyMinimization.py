import sys,os
from subprocess import PIPE, run
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from progress_bar import progress_bar

# adjust saved figure dpi
matplotlib.rcParams.update({'savefig.dpi'     : 300,
                            'xtick.labelsize' : 16,
                            'ytick.labelsize' : 16   })

# make C++ executable
def make(program):
  os.system("make PROG=VMCcomputation_{:s}.o".format(program))

# run program and extract command line output
def run_program(omega,alpha,beta,N_MC,N_burn):
  result = run("./run.x {:f} {:f} {:f} 1 {:d} {:d}".format(omega,alpha,beta,N_MC,N_burn),
               stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
  return [float(line.split("=")[-1]) for line in list(filter(None,result.stdout.split("\n")))]

# adds the correct ordinal suffix to printed numbers
def ordinal(n):  
  """
  This was shamelessly stolen with love from a nice Mr. Ben Davis on Stack Exchange.
  The function prints the correct abbreviation for ordinal numbers. Ex: 1st, 2nd, 3rd, 4th.
  Reference: https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement
  """
  return "%d%s" % (n,"tsnrhtdd"[(np.floor(n/10)%10!=1)*(n%10<4)*n%10::4])

"""
Program for performing a Variational Monte Carlo analysis for two Coulomb-interacting
electrons in a harmonic oscillator. The VMC is done using the C++ program: "VMCcomputation_Energy.cpp".

Variational wave function:

  psi1 = exp( - 0.5 * alpha * omega * ( (r_1)^2 + (r_2)^2 ) )
  
where (r_i)^2 is the norm-2 of the position vector of particle i.

The C++ programs' command line arguments:  ./run.x omega alpha beta TYPE N_MC N_burn
where:
  omega  | harmonic oscillator frequency
  alpha  | variational parameter
  beta   | variational parameter (used with psi2, see Wave2_EnergyMinimization.py)
  TYPE   | indicates whether psi1 or psi2 is used (TYPE=1)
  N_MC   | number of Monte Carlo cycles
  N_burn | number of Metropolis burn-in cycles
"""

# quick-selection
run_section1 = False  # Analyse Monte Carlo stability
run_section2 = True   # Analyse the variational parameter space w/ respect to energy minimization
run_section3 = False   # Analyse expected separation between the electrons
tag          = "2"     # extra save tag used in np.save/np.load of results
load_results = True    # if True, skip run computations

###
### Section 1: Analyse Monte Carlo stability
###

if run_section1:  
  # load results or run computations
  if load_results:
    N_MC,TrialEnergy,EnergyVariance,Acceptance = np.load("../results/wave1/stability_{:s}.npy".format(tag))
  else:
    # analysis parameters
    repetitions  = 50      # number of times the VMC is run with each N_MC value
    MinExponent  = 2       # 10^MinExponent is the smallest N_MC value
    MaxExponent  = 7       # 10^MaxExponent is the largest N_MC value
    numN_MC      = 21      # number of different N_MC values
    
    # misc VMC parameters
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
        TrialEnergy[i,j],EnergyVariance[i,j],Acceptance[i,j] = run_program(1,1,0,n_MC,0)
        sanity.update()
    N_MC = np.tile(N_MC,(repetitions,1)).T
    np.save("../results/wave1/stability_{:s}.npy".format(tag),[N_MC,TrialEnergy,EnergyVariance,Acceptance])
  
  # scatter-plot TrialEnergy, EnergyVariance and acceptance
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,8))
  
  fig.suptitle("Variational Monte Carlo Stability",fontsize=20)
  
  axes[0].scatter(N_MC,TrialEnergy,   s=2.0,c='r',marker='o')
  axes[1].scatter(N_MC,EnergyVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(N_MC,Acceptance,    s=2.0,c='r',marker='o')
  
  axes[0].set_ylabel(r"$E_T$",fontsize=18)
  axes[1].set_ylabel(r"Var$\,[E_T]$",fontsize=18)
  axes[2].set_ylabel("acceptance",fontsize=18)
  axes[2].set_xlabel("Monte Carlo cycles",fontsize=18)
  
  axes[2].set_xlim([N_MC.min()*0.9,N_MC.max()*1.1])
  
  axes[1].set_ylim([0,EnergyVariance.max()])
  axes[0].set_xscale("log")
  axes[1].set_xscale("log")
  axes[2].set_xscale("log")
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.92)
  
  fig.savefig("../results/wave1/stability_{:s}.png".format(tag))
  plt.show()


###
### Section 2: Analyse the variational parameter space w/ respect to energy minimization
###

if run_section2:
  # load results or run computations
  if load_results:
    alpha,TrialEnergy,EnergyVariance,Acceptance = np.load("../results/wave1/energy_minimization_{:s}.npy".format(tag))
  else:
    # parameters based on Section 1:
    N_burn = 100000   # Metropolis burn-in
    N_MC   = 100000   # number of Monte Carlo cycles
    
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
  
  fig.suptitle("Wave 1 Energy Minimization",fontsize=20)
  
  axes[0].scatter(alpha,TrialEnergy,   s=2.0,c='r',marker='o')
  axes[1].scatter(alpha,EnergyVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(alpha,Acceptance,    s=2.0,c='r',marker='o')
  
  axes[0].set_ylabel("trial energy",fontsize=18)
  axes[1].set_ylabel("energy variance",fontsize=18)
  axes[2].set_ylabel("acceptance",fontsize=18)
  axes[2].set_xlabel(r"$\alpha$",fontsize=18)
  
#  axes[0].set_yscale("log")
#  axes[1].set_yscale("log")
  
  axes[2].set_xlim([alpha.min()*0.9,alpha.max()*1.1])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.92)
  
  fig.savefig("../results/wave1/energy_minimization_{:s}.png".format(tag))
  plt.show()
  
  # polynomial regression of TrialEnergy
  deg          = 5
  p            = np.polyfit(alpha.ravel(),TrialEnergy.ravel(),deg)
  alpha_       = np.linspace(alpha.min()*0.95,alpha.max()*1.05,20000)
  TrialEnergy_ = np.sum([p[i]*alpha_**(len(p)-1-i) for i in range(len(p))],axis=0)
  
  print("approximate polynomial:")
  print("p(x) = ",end="",flush=True)
  [print("{:+.2f}x^{:d}".format(p[i],len(p)-1-i),end=" ",flush=True) for i in range(len(p))]
  print("")
  
  # plot polynomial fit on-top of data points
  fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
  
  ax.scatter(alpha, TrialEnergy, c="b",label="data",s=1.5)
  ax.plot(   alpha_,TrialEnergy_,c="r",label="polyfit")
  
  ax.set_title("{:s} Degree Polyfit of Wave 1 Energy Minimum".format(ordinal(deg)),fontsize=20)
  ax.set_xlabel(r"$\alpha$",fontsize=18)
  ax.set_ylabel(r"$E_T(\alpha)$",fontsize=18)
  
  fig.savefig("../results/wave1/energy_minimum_deg{:d}polyfit_{:s}.png".format(deg,tag))
  plt.show()
  
  # find energy minimum
  idx_min = np.argmin(TrialEnergy_)
  print("Section 2 Energy minimum:")
  print("E     = {:.8e}".format(TrialEnergy_[idx_min]))
  print("alpha = {:.8e}".format(alpha_[idx_min]))
  
  """
  sample run
  
  Section 2 Energy minimum:
  E     = 3.77315897e+00
  alpha = 8.79441472e-01
  """
  

###
### Section 3: Analyse expected separation between the electrons
###

if run_section3:
  # load results or run computations
  if load_results:
    omega,Separation,SeparationVariance,Acceptance = np.load("../results/wave1/expected_separation_{:s}.npy".format(tag))
  else:
    # parameters based on Sections 1 and 2:
    N_burn = 100000           # Metropolis burn-in
    N_MC   = 100000           # number of Monte Carlo cycles
    alpha  = 8.79441472e-01   # variational parameter
    
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
  
  # plot expected separation
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(10,8))
  
  fig.suptitle("Variational Monte Carlo Expected Separation Between Electrons",fontsize=20)
  
  axes[0].scatter(omega,Separation,        s=2.0,c='r',marker='o')
  axes[1].scatter(omega,SeparationVariance,s=2.0,c='r',marker='o')
  axes[2].scatter(omega,Acceptance,        s=2.0,c='r',marker='o')
  
  axes[0].set_ylabel("expected\n"+"separation",fontsize=18)
  axes[1].set_ylabel("separation\n"+"variance",fontsize=18)
  axes[2].set_ylabel("acceptance",fontsize=18)
  axes[2].set_xlabel(r"$\omega$",fontsize=18)
  
  axes[2].set_xlim([omega.min()*0.9,omega.max()*1.1])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.92)
  
  fig.savefig("../results/wave1/expected_separation_{:s}.png".format(tag))
  plt.show()
  
  
  
# clean up .o files
os.system("make clean")

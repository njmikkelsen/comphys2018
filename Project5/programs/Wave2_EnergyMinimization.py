import sys,os
from subprocess import PIPE, run
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
from progress_bar import progress_bar

# adjust saved figure dpi
matplotlib.rcParams.update({'savefig.dpi':300})

# make C++ executable
def make(program):
  os.system("make PROG=VMCcomputation_{:s}.o".format(program))

# run program and extract command line output
def run_program(omega,alpha,beta,N_MC,N_burn):
  result = run("./run.x {:f} {:f} {:f} 2 {:d} {:d}".format(omega,alpha,beta,N_MC,N_burn),
               stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
  return [float(line.split("=")[-1]) for line in list(filter(None,result.stdout.split("\n")))]

"""
Program for performing a Variational Monte Carlo analysis for two Coulomb-interacting
electrons in a harmonic oscillator. The VMC is done using the C++ programs:
"VMCcomputation_Energy.cpp" and "VMCcomputation_Separation.cpp".

Variational wave functions:

  psi1 = exp( - 0.5 * alpha * omega * ( (r_1)^2 + (r_2)^2 ) )
  psi2 = psi1 * exp( 0.5 * r12 / (1 + beta * r12) )
  
where (r_i)^2 is the norm-2 of the position vector of particle i and r12 is the separation between the particles.

The C++ programs' command line arguments:  ./run.x omega alpha beta TYPE N_MC N_burn
where:
  omega  | harmonic oscillator frequency
  alpha  | variational parameter
  beta   | variational parameter
  TYPE   | indicates whether psi1 or psi2 is used (=2)
  N_MC   | number of Monte Carlo cycles
  N_burn | number of Metropolis burn-in cycles
"""

# quick-selection
run_section1 = True    # Analyse Monte Carlo stability
run_section2 = False   # Analyse the variational parameter space w/ respect to energy
run_section3 = False   # Analyse expected separation between the electrons

###
### Section 1: Analyse Monte Carlo stability
###

if run_section1:  
  # data set parameters
  tag          = "2"    # extra save tag used in np.save/np.load of results
  load_results = True   # if True, skip run computations
  
  # load results or run computations
  if load_results:
    (alpha,beta),(TrialEnergy,EnergyVariance,Acceptance) = np.load("../results/wave2/energy_minimization_{:s}.npy".format(tag))
  else:
    # analysis parameters
    repetitions  = 10      # number of times the VMC is run with each pair of alpha and beta values
    alpha_center = 0.95    # optimal alpha parameter from psi1
    alpha_shift  = 0.10    # max shift in alpha from center value
    N_alpha      = 20      # number of alpha parameters
    beta_min     = 0.00    # minimum beta
    beta_max     = 0.80    # maximum beta
    N_beta       = 100     # number of beta parameters
    
    # misc parameters
    omega  = 1.0     # harmonic oscillator frequency
    N_burn = 10000   # Metropolis burn-in
    N_MC   = 10000   # number of Monte Carlo cycles
    
    # linearly spaced arrays of alpha and beta values
    alpha = alpha_center + np.linspace(-alpha_shift,alpha_shift,N_alpha)
    beta  = np.linspace(beta_min,beta_max,N_beta)
    
    # result arrays
    TrialEnergy    = np.zeros((N_alpha,N_beta,repetitions))
    EnergyVariance = np.zeros((N_alpha,N_beta,repetitions))
    Acceptance     = np.zeros((N_alpha,N_beta,repetitions))
    
    # run computations
    make("Energy")
    sanity = progress_bar(N_alpha*N_beta*repetitions)
    for k in range(repetitions):
      for i,a in enumerate(alpha):
          for j,b in enumerate(beta):
            TrialEnergy[i,j,k],EnergyVariance[i,j,k],Acceptance[i,j,k] = run_program(omega,a,b,N_MC,N_burn)
            sanity.update()
    TrialEnergy    = np.mean(TrialEnergy,   axis=-1)
    EnergyVariance = np.mean(EnergyVariance,axis=-1)
    Acceptance     = np.mean(Acceptance,    axis=-1)
    np.save("../results/wave2/energy_minimization_{:s}.npy".format(tag),[[alpha,beta],[TrialEnergy,EnergyVariance,Acceptance]])
  
  
  # image-plot Average Separation and Average SeparationVariance
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(5,10))
  
  fig.suptitle("Variational Monte Carlo Energy Minimization")
  
  im1 = axes[0].imshow(TrialEnergy,   cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  im2 = axes[1].imshow(EnergyVariance,cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  im3 = axes[2].imshow(Acceptance,    cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  
  axes[0].set_title("trial energy")
  axes[1].set_title("energy variance")
  axes[2].set_title("acceptance")
  
  xticks = np.linspace(0,len(beta)-1,5)
  yticks = np.linspace(0,len(alpha)-1, 5)
  
  xticklabels = ["{:4.2f}".format(n) for n in np.linspace(beta.min(), beta.max(),5)]
  yticklabels = ["{:4.2f}".format(n) for n in np.linspace(alpha.min(),alpha.max(), 5)]
  
  for i in range(3):
    axes[i].set_ylabel(r"$\alpha$")
    axes[i].set_yticks(yticks)
    axes[i].set_yticklabels(yticklabels)
  
  axes[2].set_xlabel(r"$\beta$")
  axes[2].set_xticks(xticks)
  axes[2].set_xticklabels(xticklabels)
  
  fig.colorbar(im1,ax=axes[0])
  fig.colorbar(im2,ax=axes[1])
  fig.colorbar(im3,ax=axes[2])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.9)
  
  fig.savefig("../results/wave2/energy_minimization_{:s}.png".format(tag))
  plt.show()
  
  
  






















# clean up .o files
os.system("make clean")

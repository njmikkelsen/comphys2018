import sys,os
from subprocess import PIPE, run
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from progress_bar import progress_bar
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

# adjust saved figure dpi
matplotlib.rcParams.update({'savefig.dpi'     : 300,
                            'xtick.labelsize' : 16,
                            'ytick.labelsize' : 16   })

# make C++ executable
def make(program):
  os.system("make PROG=VMCcomputation_{:s}.o".format(program))

# run program and extract command line output
def run_program(omega,alpha,beta,N_MC,N_burn):
  result = run("./run.x {:f} {:f} {:f} 2 {:d} {:d}".format(omega,alpha,beta,N_MC,N_burn),
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
electrons in a harmonic oscillator. The VMC is done using the C++ program: "VMCcomputation_Separation.cpp".

Variational wave functions:

  psi1 = exp( - 0.5 * alpha * omega * ( (r_1)^2 + (r_2)^2 ) )
  psi2 = psi1 * exp( 0.5 * r12 / (1 + beta * r12) )
  
where (r_i)^2 is the norm-2 of the position vector of particle i and r12 is the separation between the particles.

The C++ programs' command line arguments:  ./run.x omega alpha beta TYPE N_MC N_burn
where:
  omega  | harmonic oscillator frequency
  alpha  | variational parameter
  beta   | variational parameter
  TYPE   | indicates whether psi1 or psi2 is used (TYPE=2)
  N_MC   | number of Monte Carlo cycles
  N_burn | number of Metropolis burn-in cycles
"""

# quick-selection
run_section1 = False   # Analyse the variational parameter space w/ respect to energy
run_section2 = True    # Analyse expected separation between the electrons
tag          = "1"     # extra save tag used in np.save/np.load of results
load_results = True    # if True, skip run computations

###
### Section 1: Analyse the variational parameter space w/ respect to energy
###

if run_section1:  
  # load results or run computations
  if load_results:
    (alpha,beta),(TrialEnergy,EnergyVariance,Acceptance) = np.load("../results/wave2/energy_minimization_{:s}.npy".format(tag))
  else:
    # analysis parameters
    repetitions  = 10      # number of times the VMC is run with each pair of alpha and beta values
    alpha_center = 0.992   # optimal alpha parameter from psi1
    alpha_shift  = 0.02    # max shift in alpha from center value
    N_alpha      = 20      # number of alpha parameters
    beta_min     = 0.24    # minimum beta
    beta_max     = 0.28    # maximum beta
    N_beta       = 20      # number of beta parameters
    
    # misc parameters
    omega  = 1.0       # harmonic oscillator frequency
    N_burn = 100000    # Metropolis burn-in
    N_MC   = 100000    # number of Monte Carlo cycles
    
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
  fig,axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(4.6,10))
  
  fig.suptitle("Wave 2 Energy Minimization",fontsize=20)
  
  im1 = axes[0].imshow(TrialEnergy,   cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  im2 = axes[1].imshow(EnergyVariance,cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  im3 = axes[2].imshow(Acceptance,    cmap="nipy_spectral",aspect=float(len(beta))/float(len(alpha)))
  
  axes[0].set_title("trial energy",fontsize=18)
  axes[1].set_title("energy variance",fontsize=18)
  axes[2].set_title("acceptance",fontsize=18)
  
  xticks = np.linspace(0,len(beta)-1,5)
  yticks = np.linspace(0,len(alpha)-1, 5)
  
  xticklabels = ["{:4.2f}".format(n) for n in np.linspace(beta.min(), beta.max(),5)]
  yticklabels = ["{:4.2f}".format(n) for n in np.linspace(alpha.min(),alpha.max(), 5)]
  
  for i in range(3):
    axes[i].set_ylabel(r"$\alpha$",fontsize=18)
    axes[i].set_yticks(yticks)
    axes[i].set_yticklabels(yticklabels)
  
  axes[2].set_xlabel(r"$\beta$",fontsize=18)
  axes[2].set_xticks(xticks)
  axes[2].set_xticklabels(xticklabels)
  
  fig.colorbar(im1,ax=axes[0])
  fig.colorbar(im2,ax=axes[1])
  fig.colorbar(im3,ax=axes[2])
  
  fig.tight_layout()
  fig.subplots_adjust(top=0.9,left=0)
  
  fig.savefig("../results/wave2/energy_minimization_{:s}.png".format(tag))
  plt.show()
  
  # two-dimensional polynomial interpolation w/ sklearn
  deg          = 5
  polyfit_size = 1000
  
  beta_  = np.linspace(beta.min(), beta.max(), polyfit_size)
  alpha_ = np.linspace(alpha.min(),alpha.max(),polyfit_size)
  
  Alpha, Beta  = np.meshgrid(alpha, beta, copy=False)
  Alpha_,Beta_ = np.meshgrid(alpha_,beta_,copy=False)
  
  TrialEnergy_ = TrialEnergy.reshape(len(beta),len(alpha))
  
  poly = PolynomialFeatures(degree=deg)
  X    = poly.fit_transform(np.c_[Alpha.flatten(), Beta.flatten()])
  X_   = poly.fit_transform(np.c_[Alpha_.flatten(),Beta_.flatten()])
  
  model = LinearRegression(fit_intercept=False)
  model.fit(X,TrialEnergy.flatten())
  
  TrialEnergy_polyfit = (X_@model.coef_[:,None]).reshape(polyfit_size,polyfit_size)
  polyfit_residual    = (X@model.coef_[ :,None]).reshape(len(beta),len(alpha)) - TrialEnergy_
  
  print("Interpolation residual = {:e}".format(np.linalg.norm(polyfit_residual)))
  
  # plot polynomial interpolation
  fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1,figsize=(5.85,9.35))
  
  fig.suptitle("{:s} Degree Polyfit of\n".format(ordinal(deg))+"Wave 2 Energy Minimum",fontsize=20)
  
  im1 = ax1.imshow(TrialEnergy_polyfit,cmap="nipy_spectral")
  im2 = ax2.imshow(polyfit_residual,   cmap="nipy_spectral",aspect=float(len(alpha))/float(len(beta)))
  
  ax1.set_title("polynomial interpolation",fontsize=18)
  ax2.set_title("interpolation residuals", fontsize=18)
  
  ax1.set_ylabel(r"$\alpha$",fontsize=18)
  
  ax2.set_ylabel(r"$\alpha$",fontsize=18)
  ax2.set_xlabel(r"$\beta$", fontsize=18)
  
  ticks_ = np.linspace(0,polyfit_size-1,5)
  
  ax1.set_xticks(ticks_)
  ax1.set_yticks(ticks_)
  ax1.set_xticklabels(["" for _ in range(5)])
  ax1.set_yticklabels(yticklabels)
  
  ax2.set_xticks(yticks)
  ax2.set_yticks(xticks)
  ax2.set_xticklabels(xticklabels)
  ax2.set_yticklabels(yticklabels)
  
  fig.colorbar(im1,ax=ax1,fraction=0.0455)
  fig.colorbar(im2,ax=ax2,fraction=0.0455)

  fig.tight_layout()
  fig.subplots_adjust(top=0.86,right=0.825)  
  
  plt.savefig("../results/wave2/energy_minimization_deg{:d}polyfit_{:s}.png".format(deg,tag))
  plt.show()
  
  # find energy minimum
  idx_min = np.unravel_index(TrialEnergy_polyfit.argmin(),TrialEnergy_polyfit.shape)
  print("Section 2 Energy minimum:")
  print("E     = {:.8e}".format(TrialEnergy_polyfit[idx_min]))
  print("alpha = {:.8e}".format(Alpha_[idx_min]))
  print("beta  = {:.8e}".format(Beta_[idx_min]))
  
  """
  sample run
  
  Interpolation residual = 2.636224e-03
  Section 2 Energy minimum:
  E     = 3.73003209e+00
  alpha = 9.72600601e-01
  beta  = 2.40560561e-01
  """
  

###
### Section 2: Analyse expected separation between the electrons
###

if run_section2:  
  # load results or run computations
  if load_results:
    omega,Separation,SeparationVariance,Acceptance = np.load("../results/wave2/expected_separation_{:s}.npy".format(tag))
  else:
    # parameters based on Sections 1 and 2:
    N_burn = 100000           # Metropolis burn-in
    N_MC   = 100000           # number of Monte Carlo cycles
    alpha  = 9.72600601e-01   # variational parameter
    beta   = 2.40560561e-01   # variational parameter
    
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
    np.save("../results/wave2/expected_separation_{:s}.npy".format(tag),[omega,Separation,SeparationVariance,Acceptance])
  
  
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
  
  fig.savefig("../results/wave2/expected_separation_{:s}.png".format(tag))
  plt.show()

# clean up .o files
os.system("make clean")

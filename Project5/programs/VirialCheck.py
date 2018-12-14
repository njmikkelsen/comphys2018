import sys,os
from subprocess import PIPE, run
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from progress_bar import progress_bar

# adjust saved figure dpi
matplotlib.rcParams.update({'savefig.dpi'     : 300,
                            'xtick.labelsize' : 12,
                            'ytick.labelsize' : 12   })

# run program and extract command line output
def run_program(interaction,omega,alpha,beta,wave,N_MC,N_burn):
  result = run("./run.x {:d} {:f} {:f} {:f} {:d} {:d} {:d}".format(interaction,omega,alpha,beta,wave,N_MC,N_burn),
               stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
  return [float(line.split("=")[-1]) for line in list(filter(None,result.stdout.split("\n")))]


"""
Program for performing a Variational Monte Carlo analysis for two Coulomb-interacting
electrons in a harmonic oscillator. The VMC is done using the C++ program: "VMCcomputation_Virial.cpp".

Variational wave functions:

  psi1 = exp( - 0.5 * alpha * omega * ( (r_1)^2 + (r_2)^2 ) )
  psi2 = psi1 * exp( 0.5 * r12 / (1 + beta * r12) )
  
where (r_i)^2 is the norm-2 of the position vector of particle i and r12 is the separation between the particles.

The C++ programs' command line arguments:  ./run.x interaction omega alpha beta TYPE N_MC N_burn
where:
  interaction | indicates whether the Coulomb interaction is included
  omega       | harmonic oscillator frequency
  alpha       | variational parameter
  beta        | variational parameter
  TYPE        | indicates whether psi1 or psi2 is used (=2)
  N_MC        | number of Monte Carlo cycles
  N_burn      | number of Metropolis burn-in cycles
"""

# quick-selection
tag          = "2"     # extra save tag used in np.save/np.load of results
load_results = True    # if True, skip run computations

if load_results:
  (omega, \
   KineticEnergy1_int,KineticEnergy1_non,PotentialEnergy1_int,PotentialEnergy1_non, \
   KineticEnergy2_int,KineticEnergy2_non,PotentialEnergy2_int,PotentialEnergy2_non, \
   Acceptance) = \
  np.load("../results/virial/energy_{:s}.npy".format(tag))
else:
  # analysis parameters
  N_burn      = 100000   # Metropolis burn-in
  N_MC        = 100000   # number of Monte Carlo cycles
  omega_min   = 0.05     # minimum alpha
  omega_max   = 1.00     # maximum alpha
  N_omega     = 100      # number of alpha parameters
  repetitions = 5        # number of times the VMC is run with each omega value

  # wave 1 variational parameters
  alpha1 = 8.79441472e-01
  beta1  = 0.00000000e+00

  # wave 2 variational parameters
  alpha2 = 9.72600601e-01
  beta2  = 2.40560561e-01

  # logarithmically spaced array of alpha values
  omega = np.linspace(omega_min,omega_max,N_omega)

  # build C++ executable & clean up
  os.system("make PROG=VMCcomputation_Virial.o")
  os.system("make clean")
  
  # result arrays
  KineticEnergy1_int   = np.zeros((N_omega,repetitions))
  KineticEnergy1_non   = np.zeros((N_omega,repetitions))
  KineticEnergy2_int   = np.zeros((N_omega,repetitions))
  KineticEnergy2_non   = np.zeros((N_omega,repetitions))
  PotentialEnergy1_int = np.zeros((N_omega,repetitions))
  PotentialEnergy1_non = np.zeros((N_omega,repetitions))
  PotentialEnergy2_int = np.zeros((N_omega,repetitions))
  PotentialEnergy2_non = np.zeros((N_omega,repetitions))
  Acceptance           = np.zeros((N_omega,repetitions))

  # run computations
  sanity = progress_bar(4*N_omega*repetitions)
  for i,om in enumerate(omega):
    for j in range(repetitions):
      # wave 1
      KineticEnergy1_non[i,j],PotentialEnergy1_non[i,j],Acceptance[i,j] = run_program(0,om,alpha1,beta1,1,N_MC,N_burn)
      sanity.update()
      KineticEnergy1_int[i,j],PotentialEnergy1_int[i,j],Acceptance[i,j] = run_program(1,om,alpha1,beta1,1,N_MC,N_burn)
      sanity.update()
      # wave 2
      KineticEnergy2_non[i,j],PotentialEnergy2_non[i,j],Acceptance[i,j] = run_program(0,om,alpha2,beta2,2,N_MC,N_burn)
      sanity.update()
      KineticEnergy2_int[i,j],PotentialEnergy2_int[i,j],Acceptance[i,j] = run_program(1,om,alpha2,beta2,2,  N_MC,N_burn)
      sanity.update()
  omega = np.tile(omega,(repetitions,1)).T
  np.save("../results/virial/energy_{:s}.npy".format(tag), \
  [omega, \
  KineticEnergy1_int,KineticEnergy1_non,PotentialEnergy1_int,PotentialEnergy1_non, \
  KineticEnergy2_int,KineticEnergy2_non,PotentialEnergy2_int,PotentialEnergy2_non, \
  Acceptance])


# prepare <V>/<K> arrays
VK1_non = PotentialEnergy1_non/KineticEnergy1_non
VK1_int = PotentialEnergy1_int/KineticEnergy1_int
VK2_non = PotentialEnergy2_non/KineticEnergy2_non
VK2_int = PotentialEnergy2_int/KineticEnergy2_int

# plot arrays
fig = plt.figure(figsize=(5.8,4.5))
ax  = fig.add_subplot(111)

ax.set_title(r"Energy balance for 2 electrons in a harmonic oscillator",fontsize=14,pad=10)

ax.scatter(omega,VK1_non,s=1.7,c='r',label="non-int wave 1")
ax.scatter(omega,VK1_int,s=1.7,c='b',label="int wave 1")
ax.scatter(omega,VK2_non,s=1.7,c='g',label="non-int wave 2")
ax.scatter(omega,VK2_int,s=1.7,c='k',label="int wave 2")

ax.legend()

ax.set_xlim([0.0,1.0])
ax.set_ylim([1.0,4.0])
  
ax.set_xlabel(r"$\omega$",fontsize=12)
ax.set_ylabel(r"$\left\langle V\right\rangle\ \ /\ \ \left\langle K\right\rangle$",fontsize=12)

fig.subplots_adjust(top=0.92,right=0.93)

fig.savefig("../results/virial/energy_balance_{:s}.png".format(tag))
plt.show()

















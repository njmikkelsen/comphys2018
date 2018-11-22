import sys,os
import numpy as np
import matplotlib.pyplot as plt

"""
Program for analysing the 20x20 Ising system with respect to Burn-In time for the
Metropolis algorithm. The thermodynamic properties of interest are the following
expected values per spin:
  energy, energy^2, |net magnetisation|, net magnetisation and net magnetisation^2
The properties are determined using Monte Carlo integration w/ random states from
the Metropolis algorithm. The simulation is coded in the c++ file "experiment2.cpp".
"""

# program parameters
datafile        = "./output/experiment_data/experiment4"  # filename without extension for experiment data
file_extension  = "dat"     # file extension for experiment data file
N_BurnIn        = int(1e6)  # number of Monte Carlo cycles spent on the burn-in
N_MC            = int(1e6)  # number of Monte Carlo cycles spent on the integration
index           = 1         # extra index for savefig
run_experiments = False     # whether to perform the experiments before continuing the analysis
plot_results    = True      # whether to plot the results

# different tempeatures to use in simulation
temperature  = np.linspace(2.25,2.27,6)
lattice_size = [100]

temperature_str = " ".join([str(k) for k in temperature])

# compile c++ file
if run_experiments:
  os.system("g++ -std=c++14 -o run.x experiment4_nonparallel.cpp -O3")

# prepare arrays
eps  = np.zeros((len(lattice_size),len(temperature)))
absM = np.zeros((len(lattice_size),len(temperature)))
C_V  = np.zeros((len(lattice_size),len(temperature)))
chi  = np.zeros((len(lattice_size),len(temperature)))
time = np.zeros((len(lattice_size),len(temperature)))

# experiments
for i,L in enumerate(lattice_size):
  path = datafile + "_nonparallel_{:d}_{:d}_{:d}.{:s}".format(N_BurnIn,N_MC,L,file_extension)
  
  # run experiment
  if run_experiments:
    os.system("./run.x {:s} {:d} {:d} {:d} {:s}".format(path,L,N_BurnIn,N_MC,temperature_str))
  
  # load data
  eps[i],absM[i],C_V[i],chi[i],time[i] = np.loadtxt(path,usecols=(0,2,5,6,7),unpack=True)

# plot data
if plot_results:
  figs = [plt.figure() for _ in range(4)]
  axes = [fig.add_subplot(111) for fig in figs]
  
  cm     = plt.get_cmap('gist_rainbow')
  colors = [cm(1.*i/len(lattice_size)) for i in range(len(lattice_size))]
  
  for i,L in enumerate(lattice_size):
    axes[0].plot(temperature,eps[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[1].plot(temperature,absM[i],color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[2].plot(temperature,C_V[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[3].plot(temperature,chi[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
  
  for i in range(len(axes)):
    axes[i].legend(loc='best')
    axes[i].set_xlabel("temperature")
  
  axes[0].set_title("2-Dimensional Ising Model - Expected Energy")
  axes[1].set_title("2-Dimensional Ising Model - Expected Magnetisation")
  axes[2].set_title("2-Dimensional Ising Model - Heat capacity")
  axes[3].set_title("2-Dimensional Ising Model - Magnetic Susceptibility")
  
  axes[0].set_ylabel(r"$\langle\epsilon\rangle$")
  axes[1].set_ylabel(r"$\langle|M|\rangle$")
  axes[2].set_ylabel(r"$C_V$")
  axes[3].set_ylabel(r"$\chi$")
  
  figs[0].savefig("./output/figures/experiment4/energy.png")
  figs[1].savefig("./output/figures/experiment4/magnetisation.png")
  figs[2].savefig("./output/figures/experiment4/heat_capacity.png")
  figs[3].savefig("./output/figures/experiment4/magnetic_susceptibility.png")
  
  plt.show()




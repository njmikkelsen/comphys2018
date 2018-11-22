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
run_experiments = False     # whether to perform the experiments before continuing the analysis
plot_results    = False     # whether to plot the results

# experiment parameters
T_min  = 2.25     # min temperature
T_max  = 2.28     # max temperature
N_T    = 16       # number of temperature values
L      = [100]     # lattice sizes

# compile c++ file
if run_experiments:
  os.system("mpic++ -std=c++14 -o run.x experiment4.cpp -O3")

# prepare arrays
T    = np.linspace(T_min,T_max,N_T)
eps  = np.zeros((len(L),N_T))
absM = np.zeros((len(L),N_T))
C_V  = np.zeros((len(L),N_T))
chi  = np.zeros((len(L),N_T))
time = np.zeros((len(L),N_T))

# experiments
for i,l in enumerate(L):
  path = datafile + "_{:d}_{:d}_{:d}.{:s}".format(N_BurnIn,N_MC,l,file_extension)
  
  # run experiment
  if run_experiments:
    if os.path.isfile(path): os.remove(path)
    os.system("mpirun -n 2 ./run.x '{:s}' {:d} {:d} {:d} {:f} {:f}  {:d}".format(path,l,N_BurnIn,N_MC,T_min,T_max,N_T))
  
  # load data
  t,eps[i],absM[i],C_V[i],chi[i] = np.loadtxt(path,usecols=(0,1,3,6,7),unpack=True)
  
  # sort data
  argidx = np.argsort(t)
  eps[i],absM[i],C_V[i],chi[i] = eps[i,argidx],absM[i,argidx],C_V[i,argidx],chi[i,argidx]

# plot data
if plot_results:
  figs = [plt.figure() for _ in range(4)]
  axes = [fig.add_subplot(111) for fig in figs]
  
  cm     = plt.get_cmap('gist_rainbow')
  colors = [cm(1.*i/len(L)) for i in range(len(L))]
  
  for i,L in enumerate(L):
    axes[0].plot(T,eps[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[1].plot(T,absM[i],color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[2].plot(T,C_V[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
    axes[3].plot(T,chi[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(L))
  
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




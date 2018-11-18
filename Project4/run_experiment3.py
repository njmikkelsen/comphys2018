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
datafile        = "./output/experiment_data/experiment3"  # filename without extension for experiment data
file_extension  = "dat"     # file extension for experiment data file
N_BurnIn        = int(1e5)  # number of Monte Carlo cycles spent on the burn-in
N_MC            = int(1e6)  # number of Monte Carlo cycles spent on the integration
index           = 1         # extra index for savefig
run_experiments = False     # whether to perform the experiments before continuing the analysis
plot_results    = False     # whether to plot the results

# different tempeatures to use in simulation
temperature = [1.0,1.7,2.4]

# compile c++ file
if run_experiments:
  os.system("g++ -std=c++14 -o run.x experiment3.cpp -O3")

# prepare arrays
E     = np.zeros((len(temperature),N_MC+1))
E2    = np.zeros((len(temperature),N_MC+1))
eps   = np.zeros((len(temperature),N_MC+1))
P_eps = []
N     = np.linspace(1,N_MC,N_MC).astype(np.float_)

# experiments
for i,T in enumerate(temperature):
  path = datafile + "_{:d}_{:d}_{:f}.{:s}".format(N_BurnIn,N_MC,T,file_extension)
  
  # run experiment
  if run_experiments:
    os.system("./run.x {:s} {:f} {:d} {:d}".format(path,T,N_BurnIn,N_MC))

  # load data
  eps[i] = np.loadtxt(path,unpack=True)
  
  # extract energy frequency-distribution
  P_eps.append(np.unique(eps[i,1:],return_counts=True))

# plot data
if plot_results:
  fig1,ax1 = plt.subplots(nrows=3,ncols=1,sharex=True)
  fig2,ax2 = plt.subplots(nrows=3,ncols=1)

  cm     = plt.get_cmap('gist_rainbow')
  colors = [cm(1.*i/len(temperature)) for i in range(len(temperature))]

  fig1.suptitle("Evolution of Ising system energy",x=0.55,y=0.95,fontsize=20)
  fig2.suptitle("Freqency-distribution of system energy",x=0.55,y=0.95,fontsize=20)

  for i,e,P,T in zip(range(len(temperature)),eps,P_eps,temperature):
    ax1[i].plot(N,e[1:],color=colors[i],linewidth=0.4,label="T = {:.1f}".format(T))
    ax1[i].set_title("T = {:.1f}".format(T))
    ax1[i].set_ylabel(r"$\epsilon$")

    ax2[i].scatter(P[0],P[1],0.7,c="b",marker='s')
    ax2[i].set_title("T = {:.1f}".format(T))
    ax2[i].set_xlabel(r"$\epsilon$")
    ax2[i].set_ylabel("frequency")

  ax1[2].set_xlabel("number of Monte Carlo cycles")

  fig1.tight_layout()
  fig1.subplots_adjust(top=0.83)

  fig2.tight_layout()
  fig2.subplots_adjust(top=0.83)

  fig1.savefig("./output/figures/experiment3/system_energy.png")
  fig2.savefig("./output/figures/experiment3/energy_distribution.png")
  plt.show()

# print statistics
print("standard deviation in the system energy:")
for e,T in zip(eps,temperature):
  print("  T={:.1f} | std = {:12.8e}".format(T,np.std(e)))



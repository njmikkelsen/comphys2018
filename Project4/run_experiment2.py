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
datafile_noext  = "./output/experiment_data/experiment2"  # filename without extension for experiment data
file_extension  = ".dat"    # extension for datafile
N_MC            = int(1e6)  # number of Monte Carlo cycles
run_experiments = True      # whether to perform the experiments before plotting results

# different tempeatures to use in simulation
temperature = [1.0,1.7,2.4]

# compile c++ file
if run_experiments:
  os.system("g++ -std=c++14 -o run.x experiment2.cpp -O3")

# prepare arrays
E    = np.zeros((2,len(temperature),N_MC+1))
E2   = np.zeros((2,len(temperature),N_MC+1))
absM = np.zeros((2,len(temperature),N_MC+1))
M    = np.zeros((2,len(temperature),N_MC+1))
M2   = np.zeros((2,len(temperature),N_MC+1))
A    = np.zeros((2,len(temperature),N_MC+1))
N    = np.linspace(1,N_MC,N_MC).astype(np.float_)

# experiments
for i,T in enumerate(temperature):
  # run experiment
  if run_experiments:
    path_order    = datafile_noext + "_{:d}_{:f}".format(N_MC,T) + "_order"    + file_extension
    path_disorder = datafile_noext + "_{:d}_{:f}".format(N_MC,T) + "_disorder" + file_extension
    os.system("./run.x {:s} {:f} {:f} true".format( path_order,   T,N_MC))
    os.system("./run.x {:s} {:f} {:f} false".format(path_disorder,T,N_MC))

  # load data
  E[0,i],E2[0,i],absM[0,i],M[0,i],M2[0,i],A[0,i] = np.loadtxt(path_order,   unpack=True)
  E[1,i],E2[1,i],absM[1,i],M[1,i],M2[1,i],A[1,i] = np.loadtxt(path_disorder,unpack=True)


# plot data
figs = [plt.figure() for _ in range(6)]
axes = [fig.add_subplot(111) for fig in figs]

for ax in axes:
  ax.set_xlabel("number of Monte Carlo cycles")
  ax.set_xscale("log")
  #ax.set_yscale("log")

cm     = plt.get_cmap('gist_rainbow')
colors = [cm(1.*i/len(temperature)) for i in range(len(temperature))]

for i,T in enumerate(temperature):
  axes[0].plot(N,E[   0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  axes[1].plot(N,E2[  0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  axes[2].plot(N,absM[0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  axes[3].plot(N,M[   0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  axes[4].plot(N,M2[  0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  axes[5].plot(N,A[   0,i,1:],color=colors[i],linestyle='-',linewidth=0.7,label="T = {:3.1f} order".format(T))
  
  axes[0].plot(N,E[   1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))
  axes[1].plot(N,E2[  1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))
  axes[2].plot(N,absM[1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))
  axes[3].plot(N,M[   1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))
  axes[4].plot(N,M2[  1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))
  axes[5].plot(N,A[   1,i,1:],color=colors[i],linestyle='--',linewidth=0.7,label="T = {:3.1f} disorder".format(T))

for ax in axes: ax.legend(loc='best')
axes[5].set_yscale("log")

axes[0].set_title("evolution of expected energy")
axes[1].set_title("evolution of expected energy squared")
axes[2].set_title("evolution of expected magnetisation")
axes[3].set_title("evolution of expected net magnetisation")
axes[4].set_title("evolution of expected net magnetisation squared")
axes[5].set_title("evolution of the number of accepted proposals")

axes[0].set_ylabel(r"$\langle\epsilon\rangle$")
axes[1].set_ylabel(r"$\langle\epsilon^2\rangle$")
axes[2].set_ylabel(r"$\langle|M|\rangle$")
axes[3].set_ylabel(r"$\langle M\rangle$")
axes[4].set_ylabel(r"$\langle M^2\rangle$")
axes[5].set_ylabel("number of accepted proposals")

figs[0].savefig("./output/figures/experiment2/energy.png")
figs[1].savefig("./output/figures/experiment2/energy2.png")
figs[2].savefig("./output/figures/experiment2/magnetisation.png")
figs[3].savefig("./output/figures/experiment2/netmagnetisation.png")
figs[4].savefig("./output/figures/experiment2/netmagnetisation2.png")
figs[5].savefig("./output/figures/experiment2/accepted.png")

plt.show()



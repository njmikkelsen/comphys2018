import sys,os
import numpy as np
import matplotlib.pyplot as plt
from TwoByTwoSystem import compute_2by2

"""
Program for analysing the 2x2 Ising system.
The thermodynamic properties of interest are the following expected values per spin:
  energy, energy^2, |net magnetisation|, net magnetisation and net magnetisation^2
The properties are determined using Monte Carlo integration w/ random states from
the Metropolis algorithm. The simulation is coded in the c++ file "experiment1.cpp".
"""

# program parameters
datafile        = "./output/experiment_data/experiment1.dat"  # filename for experiment data
N_repeat        = 40     # number of runs for each N_MC
run_experiments = False  # whether to perform the experiments before analysis
plot_results    = False  # whether to plot the results of the experiments

# different "number of Monte Carlo cycles" to use
N_MC = [100,1000,10000,100000,1000000,10000000]

if run_experiments:
  # compile c++ file
  os.system("g++ -std=c++14 -o run.x experiment1.cpp -O3")

  # run experiments
  N_MC_str = " ".join([str(n) for n in N_MC]*N_repeat)
  os.system('./run.x {:s} {:s}'.format(datafile,N_MC_str))
  print("-"*54)
  print("All experiments finished, continuing the analysis.\n")

# load experiment data
N,E,E2,absM,M,M2,t = np.loadtxt(datafile,unpack=True)

if plot_results:
  # preparation for plots
  path       = "./output/figures/experiment1/"
  arrays     = [E,E2,absM,M,M2]
  titles     = ["energy","energy squared","magnetisation","net magnetisation","net magnetisation squared"]
  ylabels    = ["\epsilon","\epsilon^2","|M|","M","M^2"]
  savelabels = ["energy","energy2","magnetisation","netmagnetisation","netmagnetisation2"]
  
  figs = [plt.figure() for _ in range(6)]
  axes = [fig.add_subplot(111) for fig in figs]
    
  # plot expectance values
  for fig,ax,array,title,ylabel,savelabel in zip(figs,axes,arrays,titles,ylabels,savelabels):
    ax.scatter(N,array,1.4,c='r',marker='.')
    ax.set_title(r"Spread of expected {:s} per spin".format(title))
    ax.set_xlabel("number of Monte Carlo cycles")
    ax.set_ylabel(r"$\langle {:s}\rangle$".format(ylabel))
    ax.set_xscale("log")
    fig.subplots_adjust(left=0.15)
    fig.savefig(path+"{:s}.png".format(savelabel))
  
  # plot run times
  axes[5].scatter(N,t,1.4,c='r',marker='.')
  axes[5].set_title("Program runtimes for experiment 1")
  axes[5].set_xlabel("number of Monte Carlo cycles")
  axes[5].set_ylabel("time [sec]")
  axes[5].set_xscale("log")
  axes[5].set_yscale("log")
  figs[5].savefig(path+"time.png")
  
  plt.show()

# prepare for statistics
N    = [int(n) for n in N[:len(N_MC)]]
E    =    E.reshape(N_repeat,len(N_MC)).T
E2   =   E2.reshape(N_repeat,len(N_MC)).T
absM = absM.reshape(N_repeat,len(N_MC)).T
M    =    M.reshape(N_repeat,len(N_MC)).T
M2   =   M2.reshape(N_repeat,len(N_MC)).T
t    =    t.reshape(N_repeat,len(N_MC)).T

E_mean     = np.mean(E,   axis=1)
E2_mean    = np.mean(E2,  axis=1)
absM_mean  = np.mean(absM,axis=1)
M_mean     = np.mean(M,   axis=1)
M2_mean    = np.mean(M2,  axis=1)
t_mean     = np.mean(t,   axis=1)
E_error    = np.abs(np.std(E,   axis=1)/E_mean*   100)
E2_error   = np.abs(np.std(E2,  axis=1)/E2_mean*  100)
absM_error = np.abs(np.std(absM,axis=1)/absM_mean*100)
M_error    = np.abs(np.std(M,   axis=1)/M_mean*   100)
M2_error   = np.abs(np.std(M2,  axis=1)/M2_mean*  100)
t_error    = np.abs(np.std(t,   axis=1)/t_mean*   100)

E_exact,E2_exact,M_exact,absM_exact,M2_exact = compute_2by2(1.0,1.0)[:5]

# write statistics to file
with open("./output/experiment_data/experiment1_statistics.dat",'w') as File:
  File.write("OBS: This file is meant to be read without pagebreaks.\n\n")
  File.write(" N cycles |            eps            |           eps^2           |           abs(M)          |             M             |            M^2            |           t\n")
  File.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  for i in range(len(N_MC)):
    File.write( "{:9d} | {:+6.4e} +- {:8.3f} % | {:+6.4e} +- {:8.3f} % | {:+6.4e} +- {:8.3f} % | {:+6.4e} +- {:8.3f} % | {:+6.4e} +- {:8.3f} % | {:+6.4e} +- {:8.3f} %\n". \
    format(N[i],E_mean[i],E_error[i],E2_mean[i],E2_error[i],absM_mean[i],absM_error[i],M_mean[i],M_error[i],M2_mean[i],M2_error[i],t_mean[i],t_error[i]))
  File.write("  exact   | {:+6.4e}               | {:+6.4e}               | {:+6.4e}               | {:+6.4e}               | {:+6.4e}". \
    format(E_exact,E2_exact,absM_exact,M_exact,M2_exact))



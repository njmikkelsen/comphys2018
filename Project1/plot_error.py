import sys, os
import numpy as np
import matplotlib.pyplot as plt

# master parameters
plot_func     = False
plot_error    = True
plot_maxerror = True

# verify data files present in './data/'
try:
  files = sorted(os.listdir("./data/"))
  assert(len(files)>0)
except AssertionError:
  print("Error: No data files in 'data/'")
  print("Please create a data file using 'solve_equation.x'.")
  sys.exit(1)

# add data files to analysis
print("""Add a data file to the error analysis by typing in the filename
below. To continue press enter (with empty text). If you want
to analyse every file, type 'all'.
---------------------------------------------------------------
Available data files:""")
[print("  {:s}".format(f)) for f in files]
print("---------------------------------------------------------------")
added_files = []
while True:
  user_input = input(" | ").strip()
  if   user_input == "": break
  elif user_input == "all":
    del added_files; added_files = files; break
  elif user_input in files:
    added_files.append(files[files.index(user_input)])
added_files = sorted(added_files)

# count no. of data files with a given algorithm
k = [0]
for algorithm in ["LU","general","special","taylored"]:
  k.append(k[-1] + sum([algorithm in item for item in added_files]))

# verify at least one file added to analysis
if   len(added_files) == 0:
  print("Error: No data file added!")
  print("Please rerun the program and add a data file to the analysis.")
  sys.exit(1)
elif len(added_files) == len(files):
  print("Running error analysis on every available file.") 
elif len(added_files) == 1:
  print("Running error analysis using this file:")
  print(added_files[0])
else:
  print("Running error analysis using these files:")
  [print(f,end=", ") for f in added_files[:-1]]
  print(added_files[-1],".",sep="")

# read files and load all data
print("Reading data files. This may take some time...")
n, alg, T, dt, x, exact, y = [], [], [], [], [], [], []
for File,i in zip(added_files,range(len(added_files))):
  with open("./data/"+File, 'r') as F:
    data = F.readlines()
    # read initial lines and prepare arrays
    n.    append(int(  data[0].split("=")[-1].strip()))
    alg.  append(      data[1].split("=")[-1].strip())
    T.    append(float(data[2].split("=")[-1].strip()))
    dt.   append(float(data[3].split("=")[-1].strip()))
    x.    append(np.zeros(len(data[6:])))
    exact.append(np.zeros(len(data[6:])))
    y.    append(np.zeros(len(data[6:])))
    # read data
    for Line,j in zip(data[6:],range(len(data[6:]))):
      line        = Line.split("|")
      x[i][j]     = float(line[0].strip())
      exact[i][j] = float(line[1].strip())
      y[i][j]     = float(line[2].strip())
print("Done.")

# error analysis
print("Performing error analysis.")
eps, max_eps, h = [], np.zeros(len(added_files)), np.zeros(len(added_files))
for i in range(len(added_files)):
  eps.append(np.zeros(len(exact[i])))
  eps[i]     = np.abs((exact[i]-y[i])/exact[i])
  h[i]       = 1./(float(n[i])+1)
  max_eps[i] = np.max(eps[i])
print("Done.")
print("---------------------------------------------------------------")
print("Plotting various function plots and error plots.")

# verify plot at least one
if plot_func == plot_error == plot_maxerror == False:
  print("Error: All plots disabled.")
  print("Please reconfigure the master parameters.")
  sys.exit(1)

# plot solutions on top of exact
if plot_func:
  for i in range(len(added_files)):
    fig = plt.figure("Function plot", figsize=(10,8))
    ax  = fig.add_subplot(111)
    ax.set_title("Function plot\nSolution found using '{:s}' with n = {:d}".\
                 format(alg[i],n[i]), fontsize=22)
    ax.set_xlabel(r"$x$", fontsize=18); ax.set_xlim([0,1])
    ax.set_ylabel(r"$y$", fontsize=18)
    ax.plot(x[i],y[i],x[i],exact[i])
    ax.legend(["numerical","exact"])
    plt.savefig("figures/{:s}_{:d}_funcplot.png".format(alg[i],n[i]))
    plt.show()

# plot relative error evolution for each data file
if plot_error:
  for i in range(len(added_files)):
    # figure
    fig = plt.figure("Relative error", figsize=(14,8))
    # axes 1 - linear plot
    ax1 = fig.add_subplot(211)
    ax1.set_title("Relative error using '{:s}' with n = {:d}".format(alg[i],n[i]), fontsize=22)
    ax1.set_xlabel(r"$x$", fontsize=18); ax1.set_xlim([0,1])
    ax1.set_ylabel(r"$\epsilon$", fontsize=18)
    ax1.plot(x[i],eps[i])
    # axes 2 - logarithmic plot
    ax2 = fig.add_subplot(212)
    ax2.set_title("Logarithmic relative error using '{:s}' with n = {:d}".format(alg[i],n[i]), fontsize=22)
    ax2.set_xlabel(r"$x$", fontsize=18); ax2.set_xlim([0,1])
    ax2.set_ylabel(r"$\varepsilon$", fontsize=18)
    ax2.plot(x[i],eps[i]); ax2.set_yscale("log")
    # finish figure
    plt.tight_layout()
    plt.savefig("figures/{:s}_{:d}_error.png".format(alg[i],n[i]))
    plt.show()

# scatter plot logarithmic relative error against logarithmic step length
if plot_maxerror:
  for algorithm,i in zip(["LU","general","special","taylored"],range(4)):
    if (k[i+1]-k[i]) > 0:
      fig = plt.figure("Maximum relative error with respect to step length", figsize=(10,8))
      ax  = fig.add_subplot(111)
      ax.set_title("Logarithmic maximum relative error with respect to\nlogarithmic step length using '{:s}'".\
                   format(algorithm), fontsize=22)
      ax.set_xlabel(r"$\log(h)$", fontsize=18)
      ax.set_ylabel(r"$\varepsilon$", fontsize=18)
      ax.scatter(np.log10(h[k[i]:k[i+1]]),np.log10(max_eps[k[i]:k[i+1]]))
      plt.savefig("figures/{:s}_maxerror.png".format(algorithm))
      plt.show()

# print time spent
print("""
Time spent calculating the solutions:

algorithm        n        time spent
---------------------------------------""")
for i in range(len(added_files)):
  print("{:11s}|{:11d}| {:10f} sec".format(alg[i],n[i],T[i]/dt[i]))


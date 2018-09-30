import sys,os
import numpy as np
import matplotlib.pyplot as plt

# load filename
try:
  filename = sys.argv[1]
except:
  filename = input("Filename = ")

# read file if not already read
NumPyfilename = "HarmonicOscillator/NumPy/{:s}_{:s}.npy".format(*filename[:-4].split("_"))
if not os.path.exists(NumPyfilename):
  # read file
  try:
    with open("HarmonicOscillator/"+filename,'r') as File:
      data     = File.readlines()[3:]
      n        = int(  data[0].strip().split("=")[-1])
      infinity = float(data[1].strip().split("=")[-1])
      N_max    = int(  data[2].strip().split("=")[-1])
      epsilon  = float(data[3].strip().split("=")[-1])
      p        = int(  data[6].strip().split("=")[-1])
      a_max    = float(data[7].strip().split("=")[-1])
      dt       = float(data[8].strip().split("=")[-1][:-1])
      eigvals  = data[11:11+n]
      eigval   = np.array([float(k.strip()) for k in eigvals])
      eigvec = np.zeros((n,n))
      for i in range(n):
        k         = 13+n+(2+n)*i
        data_     = data[k+2:k+2+n]
        eigvec[i] = [float(k.strip()) for k in data_]
  except:
    print("Oops! Something went wrong. Does the file exist?")
    sys.exit(1)
  # save data as NumPy file
  np.save(NumPyfilename,[np.array([n,infinity,N_max,epsilon,p,a_max,dt]),eigval,eigvec])
else:
  data     = np.load(NumPyfilename)
  n        = int(  data[0][0])
  infinity = float(data[0][1])
  n_max    = int(  data[0][2])
  epsilon  = float(data[0][3])
  p        = int(  data[0][4])
  a_max    = float(data[0][5])
  dt       = float(data[0][6])
  eigval   = data[1]
  eigvec   = data[2]

# plot
x   = np.arange(0,len(eigval))
plt.figure(figsize=(14,8))
plt.scatter(x,eigval)
plt.title("Computed eigenvalues",fontsize=22)
plt.xlabel("eigenvalue index",fontsize=16)
plt.ylabel("eigenvalue",fontsize=16)
plt.savefig("HarmonicOscillator/plots/values{:d}_{:d}.png".format(n,int(infinity*10.)))
plt.show()

xi = np.linspace(0,infinity,n+2)
y  = np.zeros((4,n+2))
for i in range(4): y[i,1:-1] = eigvec[n-1-i]
plt.figure(figsize=(14,8))
for i in range(4): plt.plot(xi,y[i])
plt.title("First wavefunctions of the radial harmonic oscillator",fontsize=22)
plt.xlabel(r"$\xi$",fontsize=20)
plt.ylabel(r"$v(\xi)$",fontsize=20)
plt.legend(["n = {:d}".format(k) for k in range(4)],loc='best')
plt.savefig("HarmonicOscillator/plots/wavefunction{:d}_{:d}.png".format(n,int(infinity*10.)))
plt.show()


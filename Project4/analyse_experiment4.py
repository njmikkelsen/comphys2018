import sys,os
import numpy as np
import matplotlib.pyplot as plt

# program parameters
zoom      = True     # whether to zoom in on a temperature range
plot_poly = True     # plot polyfit results
plot_data = False    # wheter to plot data
save      = True     # whether to savefig
extra     = "zoom"   # extra label for savefig
T_min     = 2.26     # minimum tempearture in zoomed range  (also regression interval minimum)
T_max     = 2.30     # maximum temperature in zoomed range  (also regression interval maximum)
L_poly    = 60       # lattice size of plotted polyfit

# prepare lists of arrays
L    = []
T    = []
eps  = []
absM = []
C_V  = []
chi  = []

# read experiment data
for root,dirs,files in os.walk("./output/experiment_data/"):
  for file_ in files:
    if "experiment4" in file_:
      # store lattice dimension
      print(file_)
      L.append(int(file_.split("_")[-1].split(".")[0]))

      # load arrays
      data = np.loadtxt("./output/experiment_data/{:s}".format(file_),usecols=(0,1,3,6,7),unpack=True)
      
      # sort arrays
      t0     = data[0]
      argidx = np.argsort(t0)
      t      = t0[argidx]
      idx    = ((T_min <= t) & (t <= T_max)) if zoom else np.arange(len(t))
      
      # store arrays
      T.append(   data[0][argidx][idx])
      eps.append( data[1][argidx][idx])
      absM.append(data[2][argidx][idx])
      C_V.append( data[3][argidx][idx])
      chi.append( data[4][argidx][idx])
      
# make array
L = np.array(L)

# plot data
if plot_data:
  figs = [plt.figure() for _ in range(4)]
  axes = [fig.add_subplot(111) for fig in figs]

  cm     = plt.get_cmap('gist_rainbow')
  colors = [cm(1.*i/len(L)) for i in range(len(L))]

  for i,l in enumerate(L):
    axes[0].plot(T[i],eps[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(l))
    axes[1].plot(T[i],absM[i],color=colors[i],linewidth=0.7,label="L = {:d}".format(l))
    axes[2].plot(T[i],C_V[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(l))
    axes[3].plot(T[i],chi[i], color=colors[i],linewidth=0.7,label="L = {:d}".format(l))

  for i in range(len(axes)):
    axes[i].legend(loc='best')
    axes[i].set_xlabel(r"temperature [$\,k_B\,T\,/\,J\,$]")

  axes[0].set_title("2-Dimensional Ising Model - Expected Energy")
  axes[1].set_title("2-Dimensional Ising Model - Expected Magnetisation")
  axes[2].set_title("2-Dimensional Ising Model - Heat capacity")
  axes[3].set_title("2-Dimensional Ising Model - Magnetic Susceptibility")

  axes[0].set_ylabel(r"$\langle\epsilon\rangle$")
  axes[1].set_ylabel(r"$\langle|M|\rangle$")
  axes[2].set_ylabel(r"$C_V$")
  axes[3].set_ylabel(r"$\chi$")

  if save:
    figs[0].savefig("./output/figures/experiment4/energy_{:s}.png".format(extra))
    figs[1].savefig("./output/figures/experiment4/magnetisation_{:s}.png".format(extra))
    figs[2].savefig("./output/figures/experiment4/heat_capacity_{:s}.png".format(extra))
    figs[3].savefig("./output/figures/experiment4/magnetic_susceptibility_{:s}.png".format(extra))

  plt.show()


# second-order polynomial linear regression
Poly_Coeffs = [np.polyfit(T[i],C_V[i],deg=2) for i in range(len(L))]
T_C         = [-0.5*Poly_Coeffs[i][1]/Poly_Coeffs[i][0] for i in range(len(L))]  # polynomial vertices - h = -b/2a

print("Critical temperatures:")
for l,t_c in zip(L,T_C):
  print("L = {:3d} : T_C = {:f}".format(l,t_c))

# plot polyfit of L=L_poly data
if plot_poly:
  L_idx = np.argwhere(L==L_poly)[0][0]
  poly  = np.poly1d(Poly_Coeffs[L_idx])
  t     = np.linspace(T_min,T_max,1000)
  plt.plot(t,poly(t),label="Polyfit")
  plt.plot(T[L_idx],C_V[L_idx],label="Data")
  plt.title("Second-Order Polyomial Regression of Heat Capacity L = {:d}".format(L_poly))
  plt.xlabel(r"temperature [$\,k_B\,T\,/\,J\,$]")
  plt.ylabel(r"$C_V$")
  if save: plt.savefig("./output/figures/experiment4/polyfit_{:d}.png".format(L_poly))
  plt.show()

# plot T_C against inverse of L
"""
Critical temperatures found using these temperature intervals
  L =  40 : T_C = 2.300606, regression interval = [2.26,2.30]
  L =  60 : T_C = 2.292721, regression interval = [2.26,2.30]
  L =  80 : T_C = 2.290087, regression interval = [2.27,2.29]
  L = 100 : T_C = 2.283954, regression interval = [2.27,2.29]
"""

critical = np.array([2.300606,2.292721,2.290087,2.283954])
inv_L    = 1./np.sort(L)
inv_L_   = np.linspace(0,1.1*np.max(inv_L),1000)

Poly_crit = np.poly1d(np.polyfit(inv_L,critical,deg=1))

print("Final critical temperature: T = {:f}".format(Poly_crit(0)))

plt.scatter(inv_L,critical,label="data",color='red')
plt.plot(inv_L_,Poly_crit(inv_L_),"--",label="linear fit")
plt.title("Ising model critical temperature")
plt.xlabel("1/L")
plt.ylabel(r"$T_C$")
plt.legend(loc=2)
#plt.savefig("./output/figures/experiment4/critical.png")
plt.show()



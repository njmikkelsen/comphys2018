import numpy as np
import matplotlib.pyplot as plt

# parameters
J  = 1      # coupling strength in units eV
T0 = 1e-1   # initial temperature in units eV/k_B
T1 = 5e0    # final   temperature in units eV/k_B
NT = 1000   # number of temperature points

# temperature array & hyperbolic argument array
T   = np.linspace(T0,T1,NT)
ang = 8*J/T
  
# compute the partition function
Z = 12 + 4*np.cosh(ang)

# compute the energy- and net spin-variances
var_eps  = (128.*J*J/Z/Z) * (Z*np.cosh(ang) - 8*np.sinh(ang)**2)
var_spin = (32/Z) * (1 + np.exp(ang))

# compute heat capacity & magnetic susceptibility
C_V = var_eps / (T*T)
chi = var_spin / T

# plot function
fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True)

ax1.plot(T,C_V)
ax1.set_title("Heat capacity")
ax1.set_ylabel(r"$C_V$")
ax1.set_xlim([0,ax1.get_xlim()[1]])

ax2.plot(T,chi)
ax2.set_title("Magnetic Susceptibility")
ax2.set_xlabel(r"Temperature")
ax2.set_ylabel(r"$\chi$")

fig.tight_layout()
plt.subplots_adjust(top=0.85)
fig.suptitle(r"$2\times2$ Ising Model",y=0.97,x=0.54,fontsize=20)

plt.savefig("./output/figures/HeatCapacity_MagneticSusceptibility_TwoByTwo.png")
plt.show()


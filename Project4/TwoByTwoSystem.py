import numpy as np
import matplotlib.pyplot as plt

# compute quantities
def compute_2by2(temperature,J=1):
  # exponential argument array
  ang = 8*J/temperature
    
  # compute the partition function
  Z = 8 + 2*np.exp(-ang) + 2*np.exp(ang)

  # compute energy quantities
  eps     = 16*J*    (np.exp(-ang) - np.exp(ang))/Z
  eps_2   = 128*J**2*(np.exp(-ang) + np.exp(ang))/Z
  var_eps = eps_2 - eps**2

  # compute magnetic quantities
  if hasattr(Z,"__len__"):
    M = np.zeros(len(Z))
  else:
    M = 0
  absM  = (16 + 4*np.exp(ang))/Z
  M_2   = 32*(1 + np.exp(ang))/Z
  var_M = M_2 - absM**2

  # compute heat capacity and susceptibility
  C_V = var_eps/(temperature*temperature)
  chi = var_M/temperature
  
  return eps/4,eps_2/4,M/4,absM/4,M_2/4,C_V/4,chi/4


if __name__ == "__main__":
  # parameters
  J  = 1      # coupling strength in units eV
  T0 = 1e-1   # initial temperature in units eV/k_B
  T1 = 1e1    # final   temperature in units eV/k_B
  NT = 1000   # number of temperature points

  # compute
  T = np.linspace(T0,T1,NT)
  eps,eps_2,M,absM,M_2,C_V,chi = compute_2by2(T,J)

  # plot energy quantites
  fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True)

  ax1.plot(T,eps)
  ax1.set_title("Expected energy per spin")
  ax1.set_ylabel(r"$\langle\epsilon\rangle\,/\,N_{spin}$")
  ax1.set_xlim([0,ax1.get_xlim()[1]])

  ax2.plot(T,eps_2)
  ax2.set_title("Expected energy squared per spin")
  ax2.set_xlabel(r"Temperature $k_BT/J$")
  ax2.set_ylabel(r"$\langle\epsilon^2\rangle\,/\,N_{spin}$")

  fig.tight_layout()
  plt.subplots_adjust(top=0.85)
  fig.suptitle(r"$2\times2$ Ising Model",y=0.97,x=0.54,fontsize=20)

  plt.savefig("./output/figures/2by2/energy.png")
  plt.show()

  # plot magnetic quantities
  fig,(ax1,ax2,ax3) = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(8,8))

  ax1.plot(T,absM)
  ax1.set_title("Expected magnetisation per spin")
  ax1.set_ylabel(r"$\langle\|M|\rangle\,/\,N_{spin}$")
  ax1.set_xlim([0,ax1.get_xlim()[1]])

  ax2.plot(T,M)
  ax2.set_title("Expected net magnetisation per spin")
  ax2.set_xlabel(r"Temperature $k_BT/J$")
  ax2.set_ylabel(r"$\langle\ M\rangle\,/\,N_{spin}$")

  ax3.plot(T,M_2)
  ax3.set_title("Expected net magnetisation squared per spin")
  ax3.set_xlabel(r"Temperature $k_BT/J$")
  ax3.set_ylabel(r"$\langle\ M^2\rangle\,/\,N_{spin}$")

  fig.tight_layout()
  plt.subplots_adjust(top=0.85)
  fig.suptitle(r"$2\times2$ Ising Model",y=0.97,x=0.54,fontsize=20)

  plt.savefig("./output/figures/2by2/magnetisation.png")
  plt.show()

  # plot heat capacity and susceptibility function
  fig,(ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True)

  ax1.plot(T,C_V)
  ax1.set_title("Heat capacity per spin")
  ax1.set_ylabel(r"$C_V\,/\,N_{spin}$")
  ax1.set_xlim([0,ax1.get_xlim()[1]])

  ax2.plot(T,chi)
  ax2.set_title("Magnetic Susceptibility per spin")
  ax2.set_xlabel(r"Temperature $k_BT/J$")
  ax2.set_ylabel(r"$\chi\,/\,N_{spin}$")

  fig.tight_layout()
  plt.subplots_adjust(top=0.85)
  fig.suptitle(r"$2\times2$ Ising Model",y=0.97,x=0.54,fontsize=20)

  plt.savefig("./output/figures/2by2/heatcapacity_susceptibility.png")
  plt.show()



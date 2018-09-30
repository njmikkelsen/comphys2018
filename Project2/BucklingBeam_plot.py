import sys,os
import numpy as np
import matplotlib.pyplot as plt

# load filename
try:
  filename = sys.argv[1]
except:
  filename = input("Filename = ")

# read file if not already read
NumPyfilename = "BucklingBeam/NumPy/{:s}_{:s}.npy".format(*filename[:-4].split("_"))
if not os.path.exists(NumPyfilename):
  # read file
  try:
    with open("BucklingBeam/"+filename,'r') as File:
      data    = File.readlines()[3:]
      n       = int(  data[0].strip().split("=")[-1])
      N_max   = int(  data[1].strip().split("=")[-1])
      epsilon = float(data[2].strip().split("=")[-1])
      p       = int(  data[5].strip().split("=")[-1])
      a_max   = float(data[6].strip().split("=")[-1])
      dt      = float(data[7].strip().split("=")[-1][:-1])
      eigval  = np.zeros((n,3))
      for i in range(n):
        line      = data[11+i].split("|")
        eigval[i] = [float(ele.strip()) for ele in line]
      eigvec = np.zeros((n,2,n))
      for i in range(n):
        k     = 14+n+(3+n)*i
        data_ = data[k+2:k+2+n]
        for j in range(n):
          line            = data_[j].strip().split(" | ")
          eigvec[i][0][j] = float(line[0])
          eigvec[i][1][j] = float(line[1])
  except:
    print("Oops! Something went wrong. Does the file exist?")
    sys.exit(1)
  # save data as NumPy file
  np.save(NumPyfilename,[np.array([n,N_max,epsilon,p,a_max,dt]),eigval,eigvec])
else:
  data    = np.load(NumPyfilename)
  n       = int(  data[0][0])
  n_max   = int(  data[0][1])
  epsilon = float(data[0][2])
  p       = int(  data[0][3])
  a_max   = float(data[0][4])
  dt      = float(data[0][5])
  eigval  = data[1]
  eigvec  = data[2]
  
# print MSE scores
MSE_Jacobi = np.sum((eigval[:,0]-eigval[:,1])**2)
MSE_arma   = np.sum((eigval[:,0]-eigval[:,2])**2)
print("Jacobi    MSE =",MSE_Jacobi)
print("armadillo MSE =",MSE_arma)

# plot
x   = np.arange(0,len(eigval))
fig = plt.figure(figsize=(14,8))
ax1 = fig.add_subplot(211)
ax1.set_title("Relative difference between exact and computed eigenvalues",fontsize=22)
bar1 = ax1.bar(x,100*np.abs((eigval[:,0]-eigval[:,1])/eigval[:,0]))
bar2 = ax1.bar(x,100*np.abs((eigval[:,0]-eigval[:,2])/eigval[:,0]))
for i in range(len(bar1)):
  bar1[i].set_color("b")
  bar2[i].set_color("r")
ax1.set_xlabel("eigenvalue index",fontsize=16)
ax1.set_ylabel("Relative difference [%]",fontsize=16)
ax1.legend(["Jacobi","armadillo"],loc='best')
ax2 = fig.add_subplot(212)
ax2.set_title("Value of exact eigenvalues",fontsize=22)
ax2.scatter(x,eigval[:,0])
ax2.set_xlabel("eigenvalue index",fontsize=16)
ax2.set_ylabel("eigenvalue",fontsize=16)
ax2.set_xlim([0,len(eigval)])
plt.tight_layout()
plt.savefig("BucklingBeam/plots/value{:d}_{:d}.png".format(n,-int(np.log10(epsilon))))
plt.show()

xi  = np.linspace(0,1,n+2)
y1  = np.zeros((4,n+2))
y2  = np.zeros((4,n+2))
for i in range(4):
  y1[i,1:-1] = np.abs(eigvec[i,0])
  y2[i,1:-1] = np.abs(eigvec[i,1])
fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
ax1.set_title("Jacobi Eigenvectors",fontsize=22)
for i in range(4): ax1.plot(xi,y1[i])
ax1.set_xlabel(r"$\xi$",fontsize=20)
ax1.set_ylabel(r"$v(\xi)$",fontsize=20)
ax1.set_xlim([0,1])
ax2 = fig.add_subplot(212)
ax2.set_title("Armadillo Eigenvectors",fontsize=22)
for i in range(4): ax2.plot(xi,y2[i])
ax2.set_xlabel(r"$\xi$",fontsize=20)
ax2.set_ylabel(r"$v(\xi)$",fontsize=20)
ax2.set_xlim([0,1])
plt.tight_layout()
plt.savefig("BucklingBeam/plots/vector{:d}_{:d}.png".format(n,-int(np.log10(epsilon))))
plt.show()



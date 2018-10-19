import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# simple forward-Euler integration
def Euler_integration(a,t0,r0,v0,dt,N):
  t = np.linspace(t0,t0+dt*(N-1),N)
  r = np.zeros((N,*r0.shape)); r[0] = r0
  v = np.zeros((N,*v0.shape)); v[0] = v0
  for i in range(N-1):
    v[i+1] = v[i] + dt*a(r[i])
    r[i+1] = r[i] + dt*v[i]
  return t,r,v

# Velocity-Verlet integration
def VelVerlet_integration(a,t0,r0,v0,dt,N):
  t = np.linspace(t0,t0+dt*(N-1),N)
  r = np.zeros((N,*r0.shape)); r[0] = r0
  v = np.zeros((N,*v0.shape)); v[0] = v0
  for i in range(N-1):
    a_i    = a(r[i])
    r[i+1] = r[i] + dt*v[i] + 0.5*dt**2*a_i
    v[i+1] = v[i] + 0.5*dt*(a_i+a(r[i+1]))
  return t,r,v

# gravitational acceleration | AU yr^-2
def gravity(r):
  return -4*np.pi**2*r/(np.linalg.norm(r)**3)

# compare orbit plots
def compare_orbits(t,r1,r2):
  fig1 = plt.figure(figsize=(10,8))
  fig2 = plt.figure(figsize=(10,8))
  ax11 = fig1.add_subplot(211)
  ax12 = fig1.add_subplot(212)
  ax21 = fig2.add_subplot(211,projection="3d")
  ax22 = fig2.add_subplot(212,projection="3d")
  
  ax11.set_title("Planar component of orbital trajectory")
  ax11.plot(r1[:,0],r1[:,1],label="Euler")
  ax11.plot(r2[:,0],r2[:,1],label="Velocity Verlet")
  ax11.legend()
  ax11.set_xlabel("x axis")
  ax11.set_ylabel("y axis")
  
  ax12.set_title("z-component of orbital trajectory")
  ax12.plot(t,r1[:,2],label="Euler")
  ax12.plot(t,r2[:,2],label="Velocity Verlet")
  ax12.legend()
  ax12.set_xlabel("time")
  ax12.set_ylabel("z axis")
  
  ax21.set_title("Orbital trajectory - Euler")
  ax21.plot(xs=r1[:,0],ys=r1[:,1],zs=r1[:,2])
  ax21.set_xlabel("x axis")
  ax21.set_ylabel("y axis")
  ax21.set_zlabel("z axis")
  
  ax22.set_title("Orbital trajectory - Velocity Verlet")
  ax22.plot(r2[:,0],r2[:,1],r2[:,2])
  ax22.set_xlabel("x axis")
  ax22.set_ylabel("y axis")
  ax22.set_zlabel("z axis")
  
  fig1.tight_layout()
  plt.show()

# parameters
t0      = 0       # initial time          | yr
N       = 300     # number of time steps  | (N/A)
dt      = 1e-3    # time step length      | yr
phi0    = 0       # angle of init pos     | rad
R       = 1e0     # radius of rotation    | AU
speed   = 2.5e0   # initial orbital speed | AU/yr
pos_rot = True    # left/right rotation   | bool

# initial conditions
c,s = np.cos(phi0),np.sin(phi0)
r0  = R*np.array([c,s,0.2*c])
v0  = speed*np.array([s,-c,0])
if pos_rot: v0 = -v0

# integration & visualisation
t1,r1,v1 = Euler_integration(gravity,t0,r0,v0,dt,N)
t2,r2,v2 = VelVerlet_integration(gravity,t0,r0,v0,dt,N)

compare_orbits(t1,r1,r2)



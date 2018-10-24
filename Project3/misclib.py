import sys,os,shutil
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from process_data import System
try:
  import cPickle as pickle
except:
  import pickle

# function for accessing the system object for sysname
def load_system(sysname):
  with open("./results/"+sysname+"/system.pkl",'rb') as File:
    return pickle.load(File)

# function for processsing system data
def save_system(sysname):
  os.system("python3 process_data.py {:s}".format(sysname))

# build program
def build(sysname):
  path = "./results/{:s}".format(sysname)
  if os.path.exists(path):
    shutil.rmtree(path)
  os.system("make PROG={:s}.o".format(sysname))

# verify output folder
def verify_output_folder(foldername):
  if not os.path.isdir("./output"):
    os.mkdir("./output")
  if not os.path.isdir("./output/{:s}".format(foldername)):
    os.mkdir("./output/{:s}".format(foldername))

# run a simulation and extract data
def run(sysname,arguments):
  os.system("./run.x {:s}".format(arguments))
  save_system(sysname)
  return load_system(sysname)

# the orbital speed of a circular orbit about the Sun
def v_circ(R):
  return 2*np.pi/np.sqrt(R)

# escape speed of a gravitational body
def v_esc(R):
  return np.sqrt(2)*v_circ(R)

# evaluate stability measures of circular orbit
def eval_circular_stability(system,dynamic_body):
  # extract data
  x,y,z,vx,vy,vz = system.get_dynamic(dynamic_body)
  R              = np.sqrt(x**2+y**2+z**2)
  V2             = vx**2+vy**2+vz**2
  r,v            = np.array([x,y,z]).T,np.array([vx,vy,vz]).T
  r_cross_v      = np.cross(r,v)
  # compute stability measures
  std_R = np.std(R)
  std_E = np.std(0.5*V2-4*np.pi/R)
  std_l = np.std(np.linalg.norm(r_cross_v,axis=1))
  return std_R,std_E,std_l

# extract the position vector array
def extract_position(system,dynamic_body):
  x,y,z = system.get_dynamic(dynamic_body)[:3]
  return np.array([x,y,z]).T

# extract array distance from the origin
def extract_dist_to_origin(system,dynamic_body):
  r = extract_position(system,dynamic_body)
  return np.linalg.norm(r,axis=1)

def plot_2d_orbit(x,y,sysname,body):
  """
  Function for illustrating a 2-dimensional orbit given by (x,y).
  """
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  
  ax.plot(x,y)
  ax.set_title("{:s} orbit".format(body))
  ax.set_xlabel("x axis [au]")
  ax.set_ylabel("y axis [au]")
  
  fig.savefig("./output/"+sysname+"/"+body+"_2d_orbit.png")
  plt.show()


def plot_3d_orbit(t,x,y,z,sysname,body,up=3):
  """
  Function for illustrating a 3-dimensional orbit given by (x,y,z), where up=1,2,3 denotes which
  coordinate axis is considered "up". Regardless of up, the coordinates are assumed to be right-handed.
  """
  if   up == 1: X,Y,Z,labels = y,z,x,["y","z","x"]  # x is up
  elif up == 2: X,Y,Z,labels = z,x,y,["z","x","y"]  # y is up
  elif up == 3: X,Y,Z,labels = x,y,z,["x","y","z"]  # z is up (default)
  
  fig1 = plt.figure()  # decomposed orbit
  fig2 = plt.figure()  # 3D orbit
  
  ax1 = fig1.add_subplot(211)
  ax2 = fig1.add_subplot(212)
  ax3 = fig2.add_subplot(111,projection="3d")
  
  ax1.plot(X,Y)
  ax1.set_title("Planar component of {:s} orbit".format(body))
  ax1.set_xlabel(labels[0]+" axis")
  ax1.set_ylabel(labels[1]+" axis")
  
  ax2.plot(t,Z)
  ax2.set_title("Time evolution of {:s} orbit height".format(body))
  ax2.set_xlabel("time")
  ax2.set_ylabel(labels[2]+" axis")
  
  ax3.plot(X,Y,Z)
  ax3.set_title("{:s} complete orbit".format(body))
  ax3.set_xlabel(labels[0]+" axis")
  ax3.set_ylabel(labels[1]+" axis")
  ax3.set_zlabel(labels[2]+" axis")
  
  fig1.tight_layout()
  
  fig1.savefig("./results/"+sysname+"/"+body+"_3d_orbit1.png")
  fig2.savefig("./results/"+sysname+"/"+body+"_3d_orbit2.png")
  plt.show()


def plot_orbit_speed(t,vx,vy,vz,sysname,body):
  """
  Function for illustrating the orbital speed of body as a function of time.
  """
  v   = np.sqrt(vx**2+vy**2+vz**2)
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  
  ax.plot(t,v)
  ax.set_title("Time evolution of {:s} orbit speed".format(body))#,fontsize=20)
  ax.set_xlabel("time")
  ax.set_ylabel("speed")
  
  fig.savefig("./results/"+sysname+"/"+body+"_orbit_speed.png")
  plt.show()

def plot_solar_system(static,dynamic,sysname):
  """
  Function for illustrating all the orbits in a gravitational solar system.
  The z-axis is assumed to be "up".
  """
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  
  r_static,v_static = static
#  for i in range(len(r_static)):
#    ax.scatter(r_static[

"""
t,general,static,dynamic = NumPy_load("The Solar System")

x,y,z,vx,vy,vz = dynamic[0]

plot_3d_orbit(t,x,y,z,"The Solar System","Earth")
plot_orbit_speed(t,vx,vy,vz,"The Solar System","Earth")
"""


















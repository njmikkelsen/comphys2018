import sys
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
try:
  import cPickle as pickle
except:
  import pickle

# system object that is stored
class System(object):
  def __init__(self,dirpath,sysname,G,N_tot,N_dyn,N_sta,N_time,time_unit,N_int,dt):
    self.dirpath    = dirpath
    self.sysname    = sysname
    self.G          = G
    self.N_tot      = N_tot
    self.N_dyn      = N_dyn
    self.N_sta      = N_sta
    self.N_time     = N_time
    self.time_unit  = time_unit
    self.time_spent = float(N_time)/float(time_unit)
    self.N_int      = N_int
    self.dt         = dt
    self.static     = {}
    self.dynamic    = []
  def add_static(self,name,r,v):
    self.static[name] = [r,v]
  def add_dynamic(self,name):
    self.dynamic.append(name)
  def get_dynamic(self,name):
    return np.load(self.dirpath+"NumPy/{:s}.npy".format(name))

if __name__ == '__main__':
  # load directory path from the command line
  try:
    assert(" ".join(sys.argv[1:]) in os.listdir("./results"))
    dirpath = "./results/"+" ".join(sys.argv[1:])+"/"
  except:
    print("missing or invalid system name command line argument")
    sys.exit(1)

  # create system object and load data from general.dat
  with open(dirpath+"general.dat",'r') as File:
    content   = File.readlines()
    sysname   = content[0][17:-2]
    G         = float(content[1].split("=")[-1].strip())
    N_tot     = int(  content[2].split("=")[-1].strip())
    N_dyn     = int(  content[3].split("=")[-1].strip())
    N_sta     = int(  content[4].split("=")[-1].strip())
    N_time    = int(  content[5].split("=")[-1].strip())
    time_unit = float(content[6].split("=")[-1].strip())
    N_int     = int(  content[7].split("=")[-1].strip())
    dt        = float(content[8].split("=")[-1].strip())
  
  system = System(dirpath,sysname,G,N_tot,N_dyn,N_sta,N_time,time_unit,N_int,dt)

  # load the static objects from "static_bodies.dat"
  with open(dirpath+"static_bodies.dat",'r') as File:
    content = File.readlines()
    for i in range(len(content)//2):
      body = content[2*i:2*i+2]
      system.add_static(body[0][:-1],body[1].split()[0:3],body[1].split()[3:6])

  # load the dynamic objects from files of type "body.dat"
  # the data is stored in NumPy arrays, system class stores the names of these arrays
  if os.path.isdir(dirpath+"NumPy/"):
    shutil.rmtree(dirpath+"NumPy/")
  os.mkdir(dirpath+"NumPy/")
  for filename in os.listdir(dirpath):
    if not filename in ["general.dat","static_bodies.dat","NumPy"]:
      data = np.loadtxt(dirpath+"/"+filename,unpack=True)
      np.save(dirpath+"NumPy/"+filename.split(".")[0],data)
      system.add_dynamic(filename.split(".")[0])

  # save system object using pickle
  with open(dirpath+"system.pkl",'wb') as File:
    pickle.dump(system,File,-1)

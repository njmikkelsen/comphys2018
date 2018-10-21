import sys
import numpy as np

# load filepath and NumPy filename
try:
  path = " ".join(sys.argv[1:])
  with open(path,'r') as F: pass
except:
  path = input("file path = ")

directory,filename = path.split("/")[-2:]
filename           = filename.split(".")[0]

# read file
with open(path,'r') as datafile:
  data = datafile.readlines()
  r    = np.zeros((len(data),3))
  v    = np.zeros(r.shape)
  for i,dataline in zip(range(len(data)),data):
    line = [float(n) for n in dataline.split()]
    r[i] = np.array(line[0:3])
    v[i] = np.array(line[3:6])

# save data as a NumPy file
np.save("./results/"+directory+"/"+filename,[r,v])


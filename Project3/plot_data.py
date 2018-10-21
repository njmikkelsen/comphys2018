import sys, os
import numpy as np
import matplotlib.pyplot as plt

try:
  path = " ".join(sys.argv[1:])
except:
  path = input("path = ")

try:
  r,v = np.load(path)
except:
  print("Error: Could not load the data!")
  print("Have you remembered to store in NumPy arrays?")
  sys.exit(1)

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(r[:,0],r[:,1])
ax2.plot(v[:,0],v[:,1])

plt.show()


import sys
import os
import shutil

# read foldername from C++
sysname = str(sys.argv[1])

# verify the results folder exists
if not os.path.isdir("./results"):
  os.mkdir("./results")

# delete folder if it exists
if os.path.isdir("./results/"+sysname):
  shutil.rmtree("./results/"+sysname)

# create directory
os.mkdir("./results/"+sysname)

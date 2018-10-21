import sys
import os
import shutil

# read foldername from C++
sysname = str(sys.argv[1])

# delete folder if it exists
if os.path.isdir("./"+sysname):
  shutil.rmtree("./"+sysname)

# create directory
os.mkdir("./"+sysname)

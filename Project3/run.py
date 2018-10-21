import os

print("Clearing old data and running makefile...")
os.system("rm ./test.x")
os.system("make clean")
os.system("make")

print("Running executable...\n")
os.system("./test.x")

print("Plotting results...")
os.system("python3 store_data.py The Solar System/Earth.dat")
os.system("python3 plot_data.py The Solar System/Earth.npy")

import sys
import os
import shutil
import numpy
import matplotlib.pyplot as plt
from misclib import *

# read program from the command line or from input
try:
  program = sys.argv[1]
except:
  program = input("program = ")
  if not program in ["all","test_algorithms","escape_velocity","three_body_problem","solar_system"]:
    print("Error: invalid program.")
    sys.exit(1)

# clear old data
if os.path.isfile("./run.x"):
  os.system("rm ./run.x")
os.system("make clean")

# run 'test_algorithms' program
if program == "test_algorithms" or program == 'all':
  build("Earth_Sun")
  verify_output_folder("test_algorithms")
  
  # parameters
  v_circ_E = v_circ(1)  # au/yr | orbital speed of Earth in circular orbit
  t_N      = 1.         # yr    | final simulation time
  N1       = 10         # (N/A) | lower number of integration steps
  N2       = 1000       # (N/A) | upper number of integration steps
  N_sim    = 20         # (N/A) | number of simulations (for both Euler and Verlet)
  
  # setup
  N     = np.linspace(N1,N2,N_sim,endpoint=False,dtype=int)
  dt    = t_N/np.asarray(N+1,dtype=float)
  std_R = np.zeros((2,N_sim))
  std_E = np.zeros((2,N_sim))
  std_l = np.zeros((2,N_sim))
  
  # run simulations
  for i in range(N_sim):
    Forward_Euler   = run("Earth_Sun","1 Euler {:f} {:f} {:d} 2".format(v_circ_E,dt[i],N[i]))
    std_R[0,i],std_E[0,i],std_l[0,i] = eval_circular_stability(Forward_Euler,"Earth")
    
    Velocity_Verlet = run("Earth_Sun","1 Verlet {:f} {:f} {:d} 2".format(v_circ_E,dt[i],N[i]))
    std_R[1,i],std_E[1,i],std_l[1,i] = eval_circular_stability(Velocity_Verlet,"Earth")
  
  # plot results
  fig1 = plt.figure(); ax1 = fig1.add_subplot(111)  # logarithmic plot of std(R) with respect to h
  fig2 = plt.figure(); ax2 = fig2.add_subplot(111)  # logarithmic plot of std(E) with respect to h
  fig3 = plt.figure(); ax3 = fig3.add_subplot(111)  # logarithmic plot of std(l) with respect to h
  
  std            = [std_R,std_E,std_l]
  ylabels        = ["R","E","l"]
  title_labels   = ["orbital radius","energy","angular momentum"]
  savefig_labels = ["radius","energy","angmom"]
  
  for fig,ax,data,ylabel,title_label,savefig_label in zip([fig1,fig2,fig3],[ax1,ax2,ax3],std,ylabels,title_labels,savefig_labels):
    for i,s in zip([0,1],["Euler","Verlet"]):
      ax.scatter(dt,data[i],label=s)
    ax.set_title("Stability of Earth's orbit - Deviation in {:s}".format(title_label))
    ax.set_xlabel(r"time step length $\Delta t$")
    ax.set_ylabel(r"std$({{{:s}}})$".format(ylabel))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5*np.min(dt),  2*np.max(dt))
    ax.set_ylim(0.5*np.min(data),2*np.max(data))
    ax.legend(loc='best')
    fig.savefig("./output/test_algorithms/Earth_Stability_{:s}.png".format(savefig_label))
  plt.show()
  
  # plot orbit
  system = run("Earth_Sun","1 Verlet {:f} 0.001 1000 2".format(v_circ_E))
  x,y    = system.get_dynamic("Earth")[:2]

  fig    = plt.figure()
  ax     = fig.add_subplot(111)

  ax.plot(x,y)
  ax.grid(True)
  ax.set_title("Earth's circular orbit about the Sun")
  ax.set_xlabel("x [au]")
  ax.set_ylabel("y [au]")

  fig.savefig("./output/test_algorithms/Earth_orbit.png")
  plt.show()  
  
  
# run 'escape velocity' program
if program == 'escape_velocity' or program == 'all':
  build("Earth_Sun")
  verify_output_folder("escape_velocity")
  
  # parameters
  v_esc_E = v_esc(1)   # au/yr | escape speed of Earth from the Sun
  t_N     = 200.       # yr    | final simulation time
  N       = 100000     # (N/A) | number of integration steps
  N_sim   = 5          # (N/A) | number of simulations
  beta1   = 2.         # (N/A) | lower beta gravity-parameter
  beta2   = 3.         # (N/A) | upper beta gravity-parameter
  
  # setup
  dt   = t_N/float(N+1)
  beta = np.linspace(beta1,beta2,N_sim)
  time = np.linspace(0,t_N,N+1)
  R    = np.zeros((N_sim,N+1))
  
  # run simulations
  for i in range(N_sim):
    system  = run("Earth_Sun","2 Verlet {:f} {:f} {:d} {:f}".format(v_esc_E,dt,N,beta[i]))
    R[i]    = extract_dist_to_origin(system,"Earth")
  
  # plot results
  fig = plt.figure()
  ax  = fig.add_subplot(111)
  
  [ax.plot(time,r,label=r"$\beta={:3.2f}$".format(b)) for r,b in zip(R,beta)]
  ax.set_title("Escape velocity - Radial trajectory of Earth")
  ax.set_xlabel("time [yr]")
  ax.set_ylabel("distance to the Sun [au]")
  ax.legend(loc='best')
  
  fig.savefig("./output/escape_velocity/radial_trajectory.png")
  plt.show()
  

# run 'three_body_problem' program
if program == "three_body_problem" or program == 'all':
  verify_output_folder("three_body_problem")

  # parameters
  v_circ_E = v_circ(1)  # au/yr | orbital speed of Earth   in circular orbit
  t_N     = 2.          # yr    | final simulation time
  N       = 200000      # (N/A) | number of integration steps
  
  # setup
  time          = np.linspace(0,t_N,N+1)
  M_factor      = [1,10,1000]
  dt            = t_N/float(N+1)
  R             = np.zeros((4,N+1))
  Earth_orbit   = np.zeros((4,N+1,3))
  Jupiter_orbit = np.zeros((N+1,3))
  
  # run simulation without Juiter
  build("Earth_Sun")
  Two_Body_System = run("Earth_Sun","1 Verlet {:f} {:f} {:d} 2".format(v_circ_E,dt,N))
  R[0]            = extract_dist_to_origin(Two_Body_System,"Earth")
  Earth_orbit[0]  = extract_position(      Two_Body_System,"Earth")

  # run simulations with Juiter
  build("Three_Body")
  for i in range(3):    
    Three_Body_System = run("Three_Body","{:f} {:d} {:f}".format(dt,N,M_factor[i]))
    R[i+1]            = extract_dist_to_origin(Three_Body_System,"Earth")
    Earth_orbit[i+1]  = extract_position(      Three_Body_System,"Earth")
    if i==0:
      Jupiter_orbit = extract_position(Three_Body_System,"Jupiter")
  
  # plot results
  fig = plt.figure()
  ax1 = fig.add_subplot(211)
  ax2 = fig.add_subplot(212)

  [ax1.plot(time,R[i],label=r"$M_{{{:s}}}$={:d}".format("factor",M_factor[i])) for i in range(3)]
  ax1.set_title("Jupiter's perturbation of Earth's orbit- dist. to the Sun")
  ax1.set_xlabel("time [yr]")
  ax1.set_ylabel("distance to the Sun [au]")
  ax1.legend(loc='best')
  
  [ax2.plot(Earth_orbit[i,:,0],Earth_orbit[i,:,1]) for i in range(3)]
  ax2.plot(Jupiter_orbit[:,0],Jupiter_orbit[:,1],label="Jupiter")
  ax2.set_title("Jupiter's perturbation of Earth's orbit - 2d orbit")
  ax2.set_xlabel("x [au]")
  ax2.set_ylabel("y [au]")
  
  fig.tight_layout()
  fig.savefig("./output/three_body_problem/Earth_perturbed_orbit.png")
  plt.show() 
  

# run 'solar_system' program
if program == 'solar_system' or program == 'all':
    build("SolarSystem")
    verify_output_folder("solar_system")
    
    # parameters
    t_N = 1     # yr    | final simulation time
    N   = 500000 # (N/A) | number of integration steps
    
    dt = t_N/float(N+1)
    
    # run simulation
    system = run("SolarSystem","{:f} {:d}".format(dt,N))
    
    # plot orbits
    fig1 = plt.figure()  # 2d plot
    fig2 = plt.figure()  # 3d plot
    
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111,projection="3d")
    
    for i in range(system.N_tot):
      # extract orbit
      body = system.dynamic[i]
      r    = extract_position(system,body)
      
      # plot
      ax1.plot(r[:,0],r[:,1],label=body)
      if not body in ["Neptune","Uranus","Saturn"]:
        ax2.plot(r[:,0],r[:,1],r[:,2])

    ax1.set_title("2-dimensional plot of the orbits of the planets")
    ax2.set_title("3-dimensional plot of the orbits of the planets")
    
    [ax.set_xlabel("x [au]") for ax in [ax1,ax2]]
    [ax.set_ylabel("y [au]") for ax in [ax1,ax2]]
    ax2.set_zlabel("z [au]")

    ax1.legend(loc='best')
    ax1.set_xlim([-2,2])
    ax1.set_ylim([-2,2])
    
    fig1.savefig("./output/solar_system/solar_system_2d.png")
    fig2.savefig("./output/solar_system/solar_system_3d.png")
    plt.show()
    

# run 'perihelion_precession' program
if program == "perihelion_precession" or program == "all":
  build("Perihelion_Precession")
  verify_output_folder("perihelion_precession")

  # parameters
  t_N = 50.0     # yr    | final simulation time
  N   = 5000000   # (N/A) | number of integration steps

  dt = float(t_N)/float(N+1)

  # run simulation
  system = run("Mercury_Sun","{:f} {:d}".format(dt,N))
  r      = extract_position(system,"Mercury")

  
  time  = np.linspace(0,t_N,N+1)
  theta = np.arctan2(r[:,1],r[:,0])

  N_ = int(0.995*N)
  
  R = np.linalg.norm(r[N_:],axis=1)
  plt.plot(R)
  plt.show()
  
  idx = N_ + np.argmin(R)
  
  rad_to_arcsec = 206264.806
  
  print("Angular shift in Perihelion = {:f}".format(abs(theta[idx]*rad_to_arcsec)))
  
  
    
    
    
    
    
    
    
    
    
    
    
  
  

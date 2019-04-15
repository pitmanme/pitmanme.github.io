# File: youngs_modulus.py
# To calculate the Young's modulus of cylindrical objects  
# Author: Mary Pitman
# Copyright 2019 Mary Pitman in Papoian Lab, Dalal Lab
# Web site: https://pitmanme.github.io
#-------------------------------------------------------------------------------#
#               to run, type (with your file locations):                        # 
#-------------------------------------------------------------------------------#

# python min_cyl_cenpa.py 'structure file' 'trajectory file'   

# Optional arugements, must be input in order,  accepts only  first 4,5,6 etc:
# 'structure file' 'trajectory file' start end jump cutoff threshold
#-------------------------------------------------------------------------------#
import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
from pylab import *
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import align
import sys

def shifted_coord(u):
    "shift COM to origin and apply principal axes rotation matrix, R"
    shift = cap.positions - cap.center_of_mass()
    Pp= np.matmul(R, shift.T) 					 
    return Pp

#Build universe
struc = sys.argv[1]
traj = sys.argv[2]
u = MDAnalysis.Universe(struc, traj)

#input parameters
#argument order:
#	  'structure file' 'trajectory' start end jump cutoff threshold 
if len(sys.argv) > 3:
    start = int(sys.argv[3])
else:
    start = 0				#start of trajectory to analyze
if len(sys.argv) > 4:
    end = int(sys.argv[4])
else:
    end = len(u.trajectory)-1		#end of trajctory to analyze
n_frames = end - start 			#number of frames for calculation
if len(sys.argv) > 5:
    jump = int(sys.argv[5])
else:
    jump = 20				#chops n_frames into n_frames/jump for RMSF/segment 
if len(sys.argv) > 6:
    cutoff = float(sys.argv[6])
else:
    cutoff= 0.6				#cutoff in Ang for RMSF high vs low
if len(sys.argv) > 7:
    threshold = int(sys.argv[7])
else:
    threshold =10 			#number of excluded C-alpha or P based on cutoff
divisor=n_frames/jump			#how many times the tajectory is chopped 
end = jump*divisor-1 			#redefine end point to be multiple of divisor

#CENP-A residue selection, can be edited to remove regions outside tested cylinder
cap = u.select_atoms('protein and name CA') + u.select_atoms('nucleic and name P') 
n_atoms = len(cap)

#-------------------------------------------------------------------------------#
# 		       	transform coord system				        #
#-------------------------------------------------------------------------------#
#Principle axes of the moment of inertia -> transformation matrix, angles between unit vectors
pa = cap.principal_axes().T
R = pa.dot(np.identity(3))

#Iterate over all time steps, data is a 3D array [N. atoms, 3 dimensions, N. frames]
data_trans = np.array([(shifted_coord(u)) for ts in u.trajectory[start:end:jump]])
data = data_trans.T

#cylindrical coordinate transform
altered_coord = np.hypot(data[:,0,:], data[:,1,:])
r_arr = np.array([altered_coord])
r_zstack= np.swapaxes(r_arr, 0, 1)


#-------------------------------------------------------------------------------#
# 		     RMSF generation and boundary finding		        #
#-------------------------------------------------------------------------------#
#calculate RMSF and create array to append so that [x, y, z, rmsf] can be sorted in columns
rmsf_array = np.zeros((0, n_atoms))
for i in range(divisor):
   start_n= i*jump
   stop_n= (i+1)*jump
   rmsf_data = RMSF(cap, start=start_n, stop=stop_n).run()
   test= rmsf_data.rmsf
   rmsf_array = np.append(rmsf_array, [test], axis=0)
rmsf_append = np.array([rmsf_array])
rmsf_zstack= np.swapaxes(rmsf_append.T, 1, 2)

#stack the shifted coordinate data and the rmsf by appending 4th column (rmsf), depth is divisor
data_set = np.concatenate((data, rmsf_zstack), axis=1)
r_set = np.concatenate((data_set, r_zstack), axis=1)

#sort rows by z magnintude and then apply rmsf cutoffs
z_sort = np.zeros((n_atoms,4, 0))
z_min = []
z_max = []
for i in range(divisor):
    i_stack = data_set[:,:,i]
    i_2d = np.squeeze(i_stack)
    i_2d = i_2d[i_2d[:,2].argsort()]    #sorting by z
    z_cutoffs = i_2d[:,3]
    z_cutoffs[z_cutoffs > cutoff] = 1
    z_cutoffs[z_cutoffs < cutoff] = 0
    z = z_cutoffs
    for s in range(len(z)):
       if np.sum(z[0:s]) == threshold:
          z_min.append(i_2d[s,2])
          break
    z_f = z[::-1]
    for s in range(len(z)):
       if np.sum(z_f[0:s]) == threshold:
          z_max.append(i_2d[-s, 2])
          break
    test2 = np.expand_dims(i_2d, axis=2)
    z_sort = np.append(z_sort, test2, axis=2)
z_height = np.asarray(z_max) - np.asarray(z_min)
z_data = z_sort

#sort by magnitude along r axis and apply cutoff
r_sort = np.zeros((n_atoms,5, 0))
r_max = []
for i in range(divisor):
    i_stack = r_set[:,:,i]
    i_2d = np.squeeze(i_stack)
    i_2d = i_2d[i_2d[:,4].argsort()]
    r_cutoffs = i_2d[:,3]
    r_cutoffs[r_cutoffs > cutoff] = 1
    r_cutoffs[r_cutoffs < cutoff] = 0
    r = r_cutoffs
    r_f = r[::-1]
    for s in range(len(r)):
       if np.sum(r_f[0:s]) == threshold:
          r_max.append(i_2d[-s, 4])
          break
    test2 = np.expand_dims(i_2d, axis=2)
    r_sort = np.append(r_sort, test2, axis=2)
r_width = np.asarray(r_max)
r_data = r_sort

print " "
#-------------------------------------------------------------------------------#
# 			      Stats on distribution			        #
#-------------------------------------------------------------------------------#
import pandas as pd

#describe distributions for z(height), r(radius)
h = pd.Series(z_height)
e_zz = (h.std())/(h.mean())
print "cylinder height description: h_avg = %.3f Angstroms, delta_h = %.3f Angstroms" % (h.mean(), h.std())

rr = pd.Series(r_width)
e_rr = (rr.std())/(rr.mean())
print "cylinder radius description: r_avg = %.3f Angstroms, delta_r = %.3f Angstroms" % (rr.mean(), rr.std())
#-------------------------------------------------------------------------------#
#  			  Young's modulus calculation			        # 
#-------------------------------------------------------------------------------#
import math

k_bT = 300.0 * 1.38e-23
V = (rr.mean())**2 * math.pi * h.mean() * 10.0**-30
pois = 0.4
E = (k_bT * ( 1.0 - pois - 2.0*pois**2))/(V*(e_zz**2 - pois*e_zz**2 + 2.0* e_rr**2 + 4.0 * pois * e_zz * e_rr))

print " "
print "Young's modulus, E= %.3f MPa" % (E*10**-6)

#-------------------------------------------------------------------------------#
#  			   Cylinder Visualization  			        #
#-------------------------------------------------------------------------------#
from mpl_toolkits.mplot3d import Axes3D 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter graph and C-alpha and phosphate positions
X = data_set[:,0,0]
Y = data_set[:,1,0]
Z = data_set[:,2,0]-z_min[0]
ax.scatter(X, Y, Z, c='#3377ff')

# define cylinder boundaries
x=np.linspace(-1,1, 100)
z=np.linspace(0, h.mean(), 100)
Xc, Zc=np.meshgrid(x, z)
Yc = np.sqrt(1-Xc**2)

# Cylinder surface parameters
rstride = 20
cstride = 10
ax.plot_surface(51*Xc, rr.mean()*Yc, Zc, color= '#6B797C',  alpha=0.5, rstride=rstride, cstride=cstride) 
ax.plot_surface(51*Xc, -rr.mean()*Yc, Zc, color= '#6B797C',  alpha=0.3, rstride=rstride, cstride=cstride) 

#Axis labels
ax.set_xlabel("P1")
ax.set_ylabel("P2")
ax.set_zlabel("P3")

savefig('demo.pdf', transparent=True)
plt.show()

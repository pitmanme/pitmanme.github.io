from math import *
from numpy import *
from matplotlib import *
from pylab import *
import MDAnalysis
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import align
###############################################
# USAGE notes
###############################################
# You must edit lines 28 and 29 so that the variables giving a file location, give the location of the file in your directories
# Once all file locations are entered, you can run script and see graphical output, data, statistics, a few tests. 
###############################################

def shifted_coord(u):
    shift = cap.positions - cap.center_of_mass()
    Pp= np.matmul(R, shift.T) 					 #R is coordinate transform matrix
    return Pp

def unshifted_coord(u):
    shift2 = cap.positions - cap.center_of_mass()
    return shift2

#H3
traj = '/Users/mpitman/work/dt/nuc/cenpc/h3/trajectories/h3.dcd'
struc = '/Users/mpitman/work/dt/nuc/cenpc/h3/structure_files/wt_nuc_noions_wow_renumbered.gro'

#Define new homolgue universe based on the alignment
u = MDAnalysis.Universe(struc, traj)

timestep = 1
start = 0
#end = len(u.trajectory)-1
end = 400
n_frames = end - start 
jump = 20
cutoff= 0.6
threshold =10 
#n_frames = len(u.trajectory)-1						
print "n_frames: ", n_frames


#H3
cap = u.select_atoms('protein and name CA') + u.select_atoms('bynum 313:4284 and name P')  + u.select_atoms('bynum 4925:8925 and name P')

n_atoms = len(cap)

#calculate RMSF and create array to append so that [x, y, z, rmsf] can be sorted in columns
divisor=n_frames/jump

rmsf_array = np.zeros((0, n_atoms))
for i in range(divisor):
   start_n= i*jump
   stop_n= (i+1)*jump
   rmsf_data = RMSF(cap, start=start_n, stop=stop_n).run()
   test= rmsf_data.rmsf
   rmsf_array = np.append(rmsf_array, [test], axis=0)

tes =  rmsf_array
rmsf_append = np.array([tes])
rmsf_zstack= np.swapaxes(rmsf_append.T, 1, 2)
print "rmsf_array", rmsf_array
#Calculate the principle axes of the moment of inertia

before = cap.principal_axes()
p1, p2, p3 = before.T
#p1, p2, p3 = cap2.principal_axes()

#Define needed inputs for transformation matrix, all unit vectors
s1 = [1, 0, 0]
s2 = [0, 1, 0]
s3 = [0, 0, 1]
one = np.dot(p3, s1)
two = np.dot(p3, s2)
three = np.dot(p3, s3)
four = np.dot(p2, s1)
five = np.dot(p2, s2)
six = np.dot(p2, s3)
seven = np.dot(p1, s1)
eight = np.dot(p1, s2)
nine = np.dot(p1, s3)
#Define tranformation matrix, currently static
R = np.array([[one, two, three], [four, five, six], [seven, eight, nine]])
#Iterate over all time steps, data is a 3D matrix defined by (N. atoms, 3 dimensions, N. frames)
data_trans = np.array([(shifted_coord(u)) for ts in u.trajectory[start:end:jump]])
data = data_trans.T

#Iterate to get the unrotated coordinates for testing
data_orig_trans = np.array([(unshifted_coord(u)) for ts in u.trajectory[start:end:jump]])
data_orig =data_orig_trans.T
#change to cylindrical coordinates
altered_coord = np.hypot(data[:,0,:], data[:,1,:])
r_arr = np.array([altered_coord])
r_zstack= np.swapaxes(r_arr, 0, 1)
print "stack shape", r_zstack.shape

#stack the shifted coordinate data and the rmsf by appending 4th column (rmsf), depth is divisor
data_set = np.concatenate((data, rmsf_zstack), axis=1)
r_set = np.concatenate((data_set, r_zstack), axis=1)
print "r_set shape", r_set.shape

#sort by mag on z-axis
z_sort = np.zeros((n_atoms,4, 0))
z_min = []
z_max = []
for i in range(divisor):
    i_stack = data_set[:,:,i]
    i_2d = np.squeeze(i_stack)
    i_2d = i_2d[i_2d[:,2].argsort()]
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
print "z_height", z_height
z_data = z_sort

#sort by magnitude along r axis
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
print "r", r_width
r_data = r_sort
print "check", r_width[1]

#describe
import pandas as pd
h = pd.Series(z_height)
print "h", h
h_output = h.describe()
print "h_output", h_output

rr = pd.Series(r_width)
r_output = rr.describe()
print "r_output", r_output


#Preliminary plotting
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter graph
X = data_set[:,0,0]
Y = data_set[:,1,0]
Z = data_set[:,2,0]-z_min[0]
ax.scatter(X, Y, Z)


print "r",  r_width[0]
r_int = r_width[0]
z_int = z_height[0]
print "z", z_height[0]
# Cylinder
x=np.linspace(-1,1, 100)
z=np.linspace(0, z_int, 100)
Xc, Zc=np.meshgrid(x, z)
Yc = np.sqrt(1-Xc**2)

# Draw parameters
rstride = 20
cstride = 10
ax.plot_surface(51*Xc, r_int*Yc, Zc, alpha=0.3, rstride=rstride, cstride=cstride)
ax.plot_surface(51*Xc, -r_int*Yc, Zc, alpha=0.3, rstride=rstride, cstride=cstride)

ax.set_xlabel("P1")
ax.set_ylabel("P2")
ax.set_zlabel("P3")
plt.show()




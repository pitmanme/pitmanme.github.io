from math import *
from numpy import *
from matplotlib import *
from pylab import *
import MDAnalysis
import numpy
import numpy.linalg
import matplotlib.pyplot as plt
import numpy as np

import matplotlib

font = {'family' : 'arial',
    'weight' : 'medium',
    'size'   : 18}

matplotlib.rc('font', **font)
################################################################################# 

get_cmap("jet")

timestep = 1
print("start loading file ....")
#pos_ca_atoms_xyz = numpy.loadtxt('positions_ca_atoms.output', unpack=True)
pos_ca_atoms_xyz = numpy.loadtxt('p_ca_t.output', unpack=True)

print("file loaded")
length_ca_atoms_xyz = len(pos_ca_atoms_xyz)
n_frames = len(pos_ca_atoms_xyz[0,:])
print("n_frames: ", n_frames)
print("length_ca_atoms_xyz: ", length_ca_atoms_xyz)
print("pos_ca_atoms_xyz: ", shape(pos_ca_atoms_xyz))


T = pos_ca_atoms_xyz.T
print("calculating covariance matrix ...")
C = numpy.cov(T,rowvar=0)

#print C
#print shape(C)

M,L = numpy.linalg.eig(C)
print(M)
print(L)

eigenvalues = numpy.zeros(length_ca_atoms_xyz)
eig_variance = numpy.zeros(length_ca_atoms_xyz)
variance_sum = 0

for i in range(length_ca_atoms_xyz):
    eigenvalues[i] = M[i]

eigenvalues.sort()

eigenvalues = eigenvalues[::-1]

for i in range(length_ca_atoms_xyz):
    variance_sum += eigenvalues[i]
    print(i, variance_sum, eigenvalues[i])

print("total variance: ", variance_sum)

for i in range(length_ca_atoms_xyz):
    eigenvalues[i] /= variance_sum

for i in range(length_ca_atoms_xyz):
    eig_variance[i] = numpy.sum(eigenvalues[0:i+1])

print("eigenvalues (normalized): ", eigenvalues)
print("eig_variance: ", eig_variance)

eig_indices = range(length_ca_atoms_xyz)

#
# variance plot
#
plt.figure()
plt.plot(eig_indices, eigenvalues)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Principal Component')
plt.ylabel('Variance (%)')
plt.show()


n = 50


m1 = L[:,0]
m2 = L[:,1]
print('m1: ',m1)
print('m2: ',m2)

#f_pc1 = open('pc1_ca_atoms_xyz.output','w')
#f_pc2 = open('pc2_ca_atoms_xyz.output','w')
#
#for i in range(length_ca_atoms_xyz):
#    f_pc1.write(str(m1[i]) + '\n')
#    f_pc2.write(str(m2[i]) + '\n')
#
#f_pc1.close()
#f_pc2.close()

#print shape(T)
#print m1

v_0 = dot(T,m1)
v_1 = dot(T,m2)

#print v_0

kb_T = 0.0019872041*300

a = zeros([n+1,n+1])

v_0_new = array(v_0)
v_1_new = array(v_1)

v_0_new -= numpy.single(numpy.mean(v_0_new))
v_1_new -= numpy.single(numpy.mean(v_1_new))

min_0 = float(min(v_0_new))
min_1 = float(min(v_1_new))

max_0 = float(max(v_0_new))
max_1 = float(max(v_1_new))

print("PCA info [min_0, max_0, min_1, max_1] : " + "[" + str(min_0) + ", " + str(max_0) + ", " + str(min_1) + ", " + str(max_1) + "].")

dr_1 = float((max_1 - min_1)/n)
dr_0 = float((max_0 - min_0)/n)


for i in range(size(v_0_new)):
    i_pc1 = int((v_0_new[i]-min_0)/dr_0)
    i_pc2 = int((v_1_new[i]-min_1)/dr_1)
    a[i_pc1][i_pc2] += 1

v0 = linspace(min_0, max_0, n+1)
v1 = linspace(min_1, max_1, n+1) 

sum_a = a.sum()

print("sum_a: ", sum_a)

F_min = 100
F_max = 0

for i in range(n+1):
    for j in range(n+1):
        if a[i][j] != 0:
            a[i][j] = -1*kb_T*log(a[i][j]/sum_a)
            if a[i][j] < F_min:
                F_min = a[i][j]
            if a[i][j] > F_max:
                F_max = a[i][j]

for i in range(n+1):
    for j in range(n+1):
        if a[i][j] == 0:
            a[i][j] = None


print("F_min: " + str(F_min) + " kT.")
print("F_max: " + str(F_max) + " kT.")

levels = (F_max - F_min)/0.5

print("The levels are: ", levels)

a = a - F_min
neg_diff = (-1)*(F_min)
a[a==neg_diff] = None
#
# free energy projection
# vmin = 0, vmax = 3.0
fig, ax0 = plt.subplots(ncols=6, sharey=True, sharex=True)
plt.figure()
plt.contourf(v0,v1,a,20, cmap='jet', vmin=0.0, vmax=4.2)
plt.title("Test, 300 K")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.xlim((-40.0,40.0))
plt.ylim((-40.0,40.0))
plt.colorbar()
plt.show()
"""
##########################################################
#
#       Calculate basin representative snapshots
#
##########################################################
box_1_frames = []
box_1_xlimits = [0, 5]
box_1_ylimits = [8.0, 12.0]
box_2_frames = []
box_2_xlimits = [-1.0, 5.0]
box_2_ylimits = [-15.0, -12.0]
for i in range(size(v_0_new)):
    if ( (box_1_xlimits[0] < v_0_new[i] < box_1_xlimits[1]) and (box_1_ylimits[0] < v_1_new[i] < box_1_ylimits[1]) ):
        box_1_frames.append(i)
for j in range(size(v_0_new)):
    if ( (box_2_xlimits[0] < v_0_new[j] < box_2_xlimits[1]) and (box_2_ylimits[0] < v_1_new[j] < box_2_ylimits[1]) ):
        box_2_frames.append(j)
# create universe                                                                                                                                                         
gro_file = 'kac_t0.gro'
xtc_file = 'kac_trunc.xtc'
print 'hello'
u = MDAnalysis.Universe(gro_file, xtc_file)
protein_or_dna_atoms = u.selectAtoms("bynum 1:11568 or nucleic")
print protein_or_dna_atoms
# user-defined variables                                                                                                                                                  
skip = 0
cut = 0
start = cut + 1 # to account for the intial gro file structure                                                                                                            
timestep = 1 # in frames                                                                                                                                                  
N_traj = len(u.trajectory) - 1 # taking away intial gro                                                                                                                   
print "N_traj: ", N_traj # check this as trajectory length (200000 if none skipped)                                                                                        
# select C-alphas and initiate geometric average                                                                                                                          
ca_atoms = u.selectAtoms('name CA')
ca_length = len(ca_atoms)
ca_array = range(len(ca_atoms))
box1_average = np.zeros((ca_length,3))
box1_length = len(box_1_frames)
box2_average = np.zeros((ca_length,3))
box2_length = len(box_2_frames)
# calculate geometric average                                                                                                                                             
for ts in u.trajectory[start::timestep]:
    t = ts.frame - 1
    if t in box_1_frames:
        print t
        for i in ca_array:
            box1_average[i][0] += ca_atoms[i].pos[0]/float(box1_length)
            box1_average[i][1] += ca_atoms[i].pos[1]/float(box1_length)
            box1_average[i][2] += ca_atoms[i].pos[2]/float(box1_length)
for ts in u.trajectory[start::timestep]:
    t = ts.frame - 1
    if t in box_2_frames:
        print t
        for i in ca_array:
            box2_average[i][0] += ca_atoms[i].pos[0]/float(box2_length)
            box2_average[i][1] += ca_atoms[i].pos[1]/float(box2_length)
            box2_average[i][2] += ca_atoms[i].pos[2]/float(box2_length)
# calculate structure with the lowest RMSD to the geometric average                                                                                                       
rmsd_min_1 = 1000.0 # rediculously high                                                                                                                                   
rmsd_min_index_1 = 2000000 # rediculously high                                                                                                                            
rmsd_min_2 = 1000.0 # rediculously high                                                                                                                                   
rmsd_min_index_2 = 2000000 # rediculously high
for ts in u.trajectory[start::timestep]:
    t = ts.frame - 1
    if t in box_1_frames:
        print t
        rmsd = 0.0
        sum = 0.0
        for i in ca_array:
            sum += pow((box1_average[i][0] - ca_atoms[i].pos[0]), 2)
            sum += pow((box1_average[i][1] - ca_atoms[i].pos[1]), 2)
            sum += pow((box1_average[i][2] - ca_atoms[i].pos[2]), 2)
            print 'sum, ', sum
            rmsd = np.sqrt(sum)/ca_length
            if (rmsd < rmsd_min_1):
                rmsd_min_1 = rmsd
                rmsd_min_index_1 = ts.frame
for ts in u.trajectory[start::timestep]:
    t = ts.frame - 1
    if t in box_2_frames:
        print t
        rmsd = 0.0
        sum = 0.0
        for i in ca_array:
            sum += pow((box2_average[i][0] - ca_atoms[i].pos[0]), 2)
            sum += pow((box2_average[i][1] - ca_atoms[i].pos[1]), 2)
            sum += pow((box2_average[i][2] - ca_atoms[i].pos[2]), 2)
            print 'sum, ', sum
            rmsd = np.sqrt(sum)/ca_length
            if (rmsd < rmsd_min_2):
                rmsd_min_2 = rmsd
                rmsd_min_index_2 = ts.frame
# output information about the average structure and the PDB, average_structure.pdb                                                                                       
print 'box1 representative rmsd, frame_1: ', rmsd_min_1, rmsd_min_index_1
print 'box2 representative rmsd, frame_2: ', rmsd_min_2, rmsd_min_index_2
# write representative structures                                                                                                                                         
for ts in u.trajectory[start::timestep]:
    if (ts.frame == rmsd_min_index_1):
        print ts.frame
        protein_or_dna_atoms.write('box3_structure_w_dna.pdb')
    elif (ts.frame == rmsd_min_index_2):
        print ts.frame
        protein_or_dna_atoms.write('box4_structure_w_dna.pdb')
"""

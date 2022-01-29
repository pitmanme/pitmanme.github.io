from math import *
from numpy import *
from matplotlib import *
from pylab import *
import MDAnalysis
import numpy
import numpy.linalg
import matplotlib.pyplot as plt


get_cmap("jet_r")

timestep = 1
print "start loading file ...."
pos_ca_atoms_xyz = numpy.loadtxt('contr_p_ca_t.output', unpack=True)
print "file loaded"
length_ca_atoms_xyz = len(pos_ca_atoms_xyz)
n_frames = len(pos_ca_atoms_xyz[0,:])
print "n_frames: ", n_frames
print "length_ca_atoms_xyz: ", length_ca_atoms_xyz
print "pos_ca_atoms_xyz: ", shape(pos_ca_atoms_xyz)


T = pos_ca_atoms_xyz.T
print "calculating covariance matrix ..."
C = numpy.cov(T,rowvar=0)

#print C
#print shape(C)

M,L = numpy.linalg.eig(C)
print M
print L

eigenvalues = numpy.zeros(length_ca_atoms_xyz)
eig_variance_percent = numpy.zeros(length_ca_atoms_xyz)
variance_sum = 0

f_eigenvalues = open('eigenvalues.output','w')
for i in range(length_ca_atoms_xyz):
    eigenvalues[i] = M[i]
    f_eigenvalues.write(str(M[i]) + '\n')

#eigenvalues.sort()
#eigenvalues = eigenvalues[::-1]

for i in range(length_ca_atoms_xyz):
    variance_sum += eigenvalues[i]
    print i, variance_sum, eigenvalues[i]

print "total variance: ", variance_sum

#for i in range(length_ca_atoms_xyz):
#    eigenvalues[i] /= variance_sum

f_eig_percent = open('eig_percent.output','w')
for i in range(length_ca_atoms_xyz):
    eig_variance_percent[i] = numpy.sum(eigenvalues[0:i+1])/variance_sum
    f_eig_percent.write(str(eig_variance_percent[i]) + '\n')

#print "eigenvalues: ", eigenvalues
#print "eig_variance: ", eig_variance

eig_indices = range(length_ca_atoms_xyz)

"""begin
plt.figure()
plt.plot(eig_indices, eigenvalues)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Principal Component')
plt.ylabel('Variance (%)')
plt.show()
end"""

n = 50


m1 = L[:,0]
m2 = L[:,1]
m3 = L[:,2]
m4 = L[:,3]
m5 = L[:,4]
print 'm1: ',m1
print 'm2: ',m2


f_pc1 = open('pc1_ca_atoms_xyz.output','w')
f_pc2 = open('pc2_ca_atoms_xyz.output','w')
f_pc3 = open('pc3_ca_atoms_xyz.output','w')
f_pc4 = open('pc4_ca_atoms_xyz.output','w')
f_pc5 = open('pc5_ca_atoms_xyz.output','w')


for i in range(length_ca_atoms_xyz):
    f_pc1.write(str(m1[i]*sqrt(eigenvalues[0])) + '\n')
    f_pc2.write(str(m2[i]*sqrt(eigenvalues[1])) + '\n')
    f_pc3.write(str(m3[i]*sqrt(eigenvalues[2])) + '\n')
    f_pc4.write(str(m4[i]*sqrt(eigenvalues[3])) + '\n')
    f_pc5.write(str(m5[i]*sqrt(eigenvalues[4])) + '\n')

f_pc1.close()
f_pc2.close()
f_pc3.close()
f_pc4.close()
f_pc5.close()
"""begin

#print shape(T)
#print m1

v_0 = dot(T,m1)
v_1 = dot(T,m2)

#print v_0

kb_T = 0.0019872041*300

a = zeros([n+1,n+1])

v_0_new = array(v_0)
v_1_new = array(v_1)

v_0_new -= numpy.float(numpy.mean(v_0_new))
v_1_new -= numpy.float(numpy.mean(v_1_new))

min_0 = float(min(v_0_new))
min_1 = float(min(v_1_new))

max_0 = float(max(v_0_new))
max_1 = float(max(v_1_new))

print "PCA info [min_0, max_0, min_1, max_1] : " + "[" + str(min_0) + ", " + str(max_0) + ", " + str(min_1) + ", " + str(max_1) + "]." 

dr_1 = float((max_1 - min_1)/n)
dr_0 = float((max_0 - min_0)/n)



for i in range(size(v_0_new)):
    i_pc1 = int((v_0_new[i]-min_0)/dr_0)
    i_pc2 = int((v_1_new[i]-min_1)/dr_1)
    a[i_pc1][i_pc2] += 1

v0 = linspace(min_0, max_0, n+1)
v1 = linspace(min_1, max_1, n+1) 

sum_a = a.sum()

#print sum_a

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


print "F_min: " + str(F_min) + " kT."
print "F_max: " + str(F_max) + " kT."

levels = (F_max - F_min)/0.5

print levels

a = a - F_min
neg_diff = (-1)*(F_min)
a[a==neg_diff] = None

# vmin = 0, vmax = 3.0
#fig, ax0 = plt.subplots(ncols=6, sharey=True, sharex=True)
plt.figure()
plt.contourf(v0,v1,a,20, cmap='jet_r', vmin=0.0, vmax=3.5)
plt.title("CENP-A 8mer, 300 K")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.xlim((-30.0,30.0))
plt.ylim((-30.0,30.0))
plt.colorbar()
plt.show()
end"""

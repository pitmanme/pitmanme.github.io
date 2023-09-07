from math import *
from numpy import *
from matplotlib import *
from pylab import *
import MDAnalysis
import numpy
import numpy.linalg
import matplotlib.pyplot as plt
import io


from MDAnalysis.analysis import pca, align
import nglview as nv


pdb = '../md1.pdb'
trj = '../md1_backbone.xtc'

# Select atoms
u = MDAnalysis.Universe(pdb,trj)
ca_atoms = u.select_atoms('name CA')
print("len(ca_atoms): ", len(ca_atoms))
print(ca_atoms)

# Trajectory parsing details
cut = 0 # times 25(frames skipped)*2(ps/frame) = 600000 ps cut
start = cut + 1
timestep = 1
n_frames = len(u.trajectory)-1						
print("n_frames: ", n_frames)

'''
# Calculate the center of mass for each timestep
tot_com = numpy.zeros((n_frames,3))
for ts in u.trajectory[start::timestep]:
    tot_com[ts.frame-1,:] = ca_atoms.center_of_mass()
'''
'''
# calculate PCA input file, CA positions with respect to COM
f_ca_atoms = open('p_ca_t.output','w')
for ts in u.trajectory[start:15:timestep]:
    pos = ca_atoms.positions
    com = ca_atoms.center_of_mass()
    dif = pos - com
    np_2_str = numpy.array2string(dif)
    position_string = str(np_2_str).replace('[','').replace(']','')
    #position_string.replace("\n", "")
    f_ca_atoms.write(position_string)
    f_ca_atoms.write('iter')
    f_ca_atoms.write('\n')
f_ca_atoms.close()
'''
''' Makes the files for other scripts
# calculate PCA input file, CA positions with respect to COM
f_ca_atoms = open('p_ca_t.output','w')
for ts in u.trajectory[start::timestep]:
    pos = ca_atoms.positions
    com = ca_atoms.center_of_mass()
    dif = pos - com
    dif_flat = dif.ravel()
    print(dif_flat)
    for i in dif_flat:
        position_str = str(i) + ' '
        f_ca_atoms.write(position_str)
    f_ca_atoms.write('\n')
f_ca_atoms.close()
'''



# New PC Analysis playing around using mdanalysis

aligner = align.AlignTraj(u, u, select='name CA', in_memory=True).run()
'''
pc = pca.PCA(u, select='name CA',
             align=True, mean=None,
             n_components=None).run()
             
n_pcs = np.where(pc.cumulated_variance > 0.95)[0][0]
print(n_pcs)
'''



'''
timestep = 1
pos_ca_atoms_xyz = numpy.loadtxt('p_ca_t.output', unpack=True)
length_ca_atoms_xyz = len(pos_ca_atoms_xyz)
n_frames = len(pos_ca_atoms_xyz[0,:])

T = pos_ca_atoms_xyz.T
print("calculating covariance matrix ...")
C = numpy.cov(T,rowvar=0)
print(f'C.shape is {C.shape}')

#M,L = numpy.linalg.eig(C)

# First question, are the result different?

print("The eivenvalues from MDAnalysis:")
print(pc.results.variance)
'''


# other way, seems best
ca_xyz = np.array([ca_atoms.positions.copy() for ts in u.trajectory])
# fortran (row-major) memory alignment to increase cache locallity
ca_xyz = ca_xyz.reshape(u.trajectory.n_frames, ca_atoms.n_atoms * 3, order='F')

x = ca_xyz - ca_xyz.mean(0)
cov = np.cov(x, rowvar=0)
e_vals, e_vecs = np.linalg.eig(cov)
# numpy sort from smallest to largest by default. We want it the other way around
sort_idx = np.argsort(e_vals)[::-1]
variance = e_vals[sort_idx]
PCs = e_vecs[:, sort_idx]

PC_projection = np.dot(x, PCs)

#m = PCs
#n = 50
n = 25

#m1 = m[:,0]
#m2 = m[:,1]

# coordinates projected onto PCs
#v_0 = dot(x,m1)
#v_1 = dot(x,m2)
v_0 = PC_projection[:, 0]
v_1 = PC_projection[:, 1]
v_2 = PC_projection[:, 2]

#print v_0
temp = 300
kb_T = 0.0019872041*temp

# Making an array of zeros, why
a = zeros([n+1,n+1])

v_0_new = array(v_0)
v_1_new = array(v_1)
v_2_new = array(v_2)

# Why are we doing this?
v_0_new -= numpy.single(numpy.mean(v_0_new))
v_1_new -= numpy.single(numpy.mean(v_1_new))
v_2_new -= numpy.single(numpy.mean(v_2_new))

min_0 = float(min(v_0_new))
min_1 = float(min(v_1_new))
min_2 = float(min(v_2_new))

max_0 = float(max(v_0_new))
max_1 = float(max(v_1_new))
max_2 = float(max(v_2_new))

print("PCA info [min_0, max_0, min_1, max_1] : " + "[" + str(min_0) + ", " + str(max_0) + ", " + str(min_1) + ", " + str(max_1) + "].")

# We are dividing the ranges by 50. Why?
dr_1 = float((max_1 - min_1)/n)
dr_0 = float((max_0 - min_0)/n)
dr_2 = float((max_2 - min_2)/n)

# tabulating how many hits you have in a box size 1
for i in range(size(v_0_new)):
    i_pc1 = int((v_0_new[i]-min_0)/dr_0)
    i_pc2 = int((v_1_new[i]-min_1)/dr_1)
    #i_pc3 = int((v_2_new[i]-min_2)/dr_2)
    a[i_pc1][i_pc2] += 1

v0 = linspace(min_0, max_0, n+1)
v1 = linspace(min_1, max_1, n+1)
v2 = linspace(min_2, max_2, n+1)

sum_a = a.sum()

#print sum_a

F_min = 100
F_max = 0

# If there is data calculate height
for i in range(n+1):
    for j in range(n+1):
        if a[i][j] != 0:
            a[i][j] = -1*kb_T*log(a[i][j]/sum_a)
            if a[i][j] < F_min:
                F_min = a[i][j]
            if a[i][j] > F_max:
                F_max = a[i][j]

# If there is no data, make nan
for i in range(n+1):
    for j in range(n+1):
        if a[i][j] == 0:
            a[i][j] = None


print("F_min: " + str(F_min) + " kT.")
print("F_max: " + str(F_max) + " kT.")

levels = (F_max - F_min)/0.5

print("levels: ", levels)

a = a - F_min
neg_diff = (-1)*(F_min)
a[a==neg_diff] = None
#
# free energy projection
# vmin = 0, vmax = 3.0
#fig, ax0 = plt.subplots(ncols=6, sharey=True, sharex=True)

'''
#---------------------#
# Make a 3D scatter plot
#---------------------#

a3d = zeros([n+1, n+1, n+1])

for i in range(size(v_0_new)):
    i_pc1 = int((v_0_new[i]-min_0)/dr_0)
    i_pc2 = int((v_1_new[i]-min_1)/dr_1)
    i_pc3 = int((v_2_new[i]-min_2)/dr_2)
    a3d[i_pc1][i_pc2][i_pc3] += 1

sum_a3d = a3d.sum()
print(sum_a3d)

# If there is data calculate height
for i in range(n+1):
    for j in range(n+1):
        for k in range(n+1):
            if a3d[i][j][k] != 0:
                a3d[i][j][k] = -1*kb_T*log(a3d[i][j][k]/sum_a3d)
                if a3d[i][j][k] < F_min:
                    F_min = a3d[i][j][k]
                if a3d[i][j][k] > F_max:
                    F_max = a3d[i][j][k]

# If there is no data, make nan
for i in range(n+1):
    for j in range(n+1):
        for k in range(n+1):
            if a3d[i][j][k] == 0:
                a3d[i][j][k] = None

print(a3d.shape)
a3d.ravel()
print(a3d.shape)
a3d = a3d[~numpy.isnan(a3d)]
print(a3d.shape)

print("F_min: " + str(F_min) + " kT.")
print("F_max: " + str(F_max) + " kT.")


a3d = a3d - F_min
neg_diff2 = (-1)*(F_min)
a3d[a3d==neg_diff2] = None

number_of_datapoints = 1842
count_min = 1
count_max = 187
data = np.random.randint(count_min, count_max, number_of_datapoints) # these are your counts
print(data)
print(data.shape)

#%% Create Color Map
colormap = plt.get_cmap("YlOrRd")
norm = matplotlib.colors.Normalize(vmin=min(data), vmax=max(data))

#%% 3D Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(v_0, v_1, v_2, s=10, c=colormap(norm(data)), marker='.')
plt.show()


ax.set_xlabel('PC1 Label')
ax.set_ylabel('PC2 Label')
ax.set_zlabel('PC3 Label')

plt.show()
'''


#print(v0)
#print(v1)
#print("a after updating: ", a)
plt.figure()
# contourf([X, Y,] Z, [levels], **kwargs)
plt.contourf(v0,v1,a,20, cmap='jet', extend='both') #vmin=0.0, vmax=4.2,)
plt.title("Test, 300 K")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.xlim((-25.0,25.0))
plt.ylim((-25.0,25.0))
plt.colorbar()
plt.show()

# Eigenvalues plot
plt.figure()
plt.plot(sort_idx, variance)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Principal Component')
plt.ylabel('Variance (%)')
plt.show()



backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print('There are {} backbone atoms in the analysis'.format(n_bb))
#print(pc.p_components.shape)
#print("principal components", pc.p_components)
#print(pc.variance)

plt.plot(pc.cumulated_variance[:10])
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance')
plt.show()

'''
import seaborn as sns
import pandas as pd



# What does this transformed do?
transformed = pc.transform(backbone, n_components=3)
df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(3)])


df['Time (ps)'] = df.index * u.trajectory.dt
print(df.head)

g = sns.PairGrid(df, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df)))
fig = g.map(plt.scatter, marker='.')
fig.savefig("out.png")
#plt.show()


pc1 = pc.p_components[:, 0]
trans1 = transformed[:, 0]
projected = np.outer(trans1, pc1) + pc.mean.flatten()
coordinates = projected.reshape(len(trans1), -1, 3)

proj1 = MDAnalysis.Merge(backbone)
proj1.load_new(coordinates, order="fac")
view = nv.show_mdanalysis(proj1.atoms)



from nglview.contrib.movie import MovieMaker
movie = MovieMaker(view, output='pc1.gif', in_memory=True)
movie.make()
'''


#
# code to analyze how tetramer angles fluctuate with PCs
#
# David Winogradoff
# 01.16.2015
#

#import modules
from math import *
from numpy import *
#from matplotlib import *
#from pylab import *
import MDAnalysis
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import math

#define helper functions: center of mass, and angle between 3 points (both in 3D)
pc_number = 1
def center_of_mass(atoms):
    com = np.zeros(3)
    for i in range(len(atoms)):
        com[0] += atoms.position[0]
        com[1] += atoms.position[1]
        com[2] += atoms.position[2]
    com[0] /= np.float(len(atoms))
    com[1] /= np.float(len(atoms))
    com[2] /= np.float(len(atoms))
    return com


def distance_between(ax,ay,az,bx,by,bz):
    a_to_b_sum = pow(ax-bx, 2) + pow(ay-by, 2) + pow(az-bz, 2)
    a_to_b_length = math.sqrt(a_to_b_sum)
    return a_to_b_length

def angle_between_3(p1x,p1y,p1z,vx,vy,vz,p2x,p2y,p2z):
    p1_to_v_sum = pow(p1x-vx, 2) + pow(p1y-vy, 2) + pow(p1z-vz, 2)
    p1_to_v_length = math.sqrt(p1_to_v_sum)
    p2_to_v_sum = pow(p2x-vx, 2) + pow(p2y-vy, 2) + pow(p2z-vz, 2)
    p2_to_v_length = math.sqrt(p2_to_v_sum)
    p2_to_p1_sum = pow(p2x-p1x, 2) + pow(p2y-p1y, 2) + pow(p2z-p1z, 2)
    p2_to_p1_length = math.sqrt(p2_to_p1_sum)
    angle_numerator = pow(p1_to_v_length, 2) + pow(p2_to_v_length, 2) - pow(p2_to_p1_length, 2)
    angle_denominator = 2.0*(p1_to_v_length)*(p2_to_v_length)
    angle_measure = math.degrees(math.acos(angle_numerator/angle_denominator))
    return angle_measure

#define or load necessary files
pc_1_column_array = "pc" + str(pc_number) + "_ca_atoms_xyz.output"
pc_1 = np.loadtxt(pc_1_column_array, dtype=np.complex128)
gro_file = '/Users/mpitman/work/dt/nuc/cenpc/contr/structure_files/cenpc_contr_t0.gro'
#xtc_file = 'wt_8mer_wow_a1_to_m13_nojump_bb_fit_600ns_to_end.xtc'
pdb_file = 'cenpc_avg_idr_dna_renum.pdb'

#create universe, and set user-defined parameters
u = MDAnalysis.Universe(pdb_file)
skip = 0
cut = 0
start = cut + 1
timestep = 1
n_frames = len(u.trajectory) - 1
multiplier = 5.0
multiplier_num = 50 #number of values desired in multiplier array
multiplier_array_a = np.linspace(-1.0*multiplier, 1.0*multiplier, num = multiplier_num)
multiplier_array_b = np.linspace(1.0*multiplier, -1.0*multiplier, num = multiplier_num)
multiplier_list = []
for a in range(len(multiplier_array_a)):
    multiplier_list.append(multiplier_array_a[a])
for b in range(1,len(multiplier_array_b)):
    multiplier_list.append(multiplier_array_b[b])

multiplier_list = 3*multiplier_list

#------------------------------------------------------------------------#

#Need to define dimer units and dna, this is with idrs removed as in 3wtp_3an2_align.xlsx
ca_atoms = u.select_atoms('name CA or name P')
AB_ca_atoms = ca_atoms[0:157] #cenpa/h4
CD_ca_atoms = ca_atoms[157:334] #h2a/h2b
E= ca_atoms[593:672] #cenpa
F= ca_atoms[334:413] #h4'
EF_ca_atoms = E+F
GH_ca_atoms = ca_atoms[413:593] #h2a'/h2b'
IJ_ca_atoms = ca_atoms[672:864] #dna
#rep_frame = 78893 #determined by another code#rep_frame = 78893 #determined by another code



residue_array = range(1,len(ca_atoms)+1)
w = MDAnalysis.Writer("contr_nuc_pc" + str(pc_number)+ "_6_loops.xtc", u.trajectory.n_atoms)
alpha_array = []
beta1_array = []
beta2_array = []
gamma_array = []
AB_to_EF = []
AB_to_CD = []
EF_to_GH = []
CD_to_GH = []
IJ_to_AB = []
IJ_to_GH = []
for ts in u.trajectory:
    print ts.frame
    for m in range(len(multiplier_list)):
        tot_com, AB_com, CD_com, EF_com, GH_com, IJ_com = [np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3)]
        for r in residue_array:
            print r
            r_string = "resid " + str(r)
            r_atoms = u.select_atoms(r_string)
            #print list(r_atoms)
            r_length = len(r_atoms)
            for n in range(r_length):
                r_atoms[n].position[0] += multiplier_list[m]*pc_1[3*(r-1)]
                r_atoms[n].position[1] += multiplier_list[m]*pc_1[3*(r-1)+1]
                r_atoms[n].position[2] += multiplier_list[m]*pc_1[3*(r-1)+2]
        
        w.write(ts)

        for ti in range(len(ca_atoms)):
            tot_com[0] += ca_atoms[ti].position[0]/np.float(len(ca_atoms))
            tot_com[1] += ca_atoms[ti].position[1]/np.float(len(ca_atoms))
            tot_com[2] += ca_atoms[ti].position[2]/np.float(len(ca_atoms))
        
        for ai in range(len(AB_ca_atoms)):
            AB_com[0] += (AB_ca_atoms[ai].position[0])/np.float(len(AB_ca_atoms))
            AB_com[1] += (AB_ca_atoms[ai].position[1])/np.float(len(AB_ca_atoms))
            AB_com[2] += (AB_ca_atoms[ai].position[2])/np.float(len(AB_ca_atoms))
        
        for ci in range(len(CD_ca_atoms)):
            CD_com[0] += (CD_ca_atoms[ci].position[0])/np.float(len(CD_ca_atoms))
            CD_com[1] += (CD_ca_atoms[ci].position[1])/np.float(len(CD_ca_atoms))
            CD_com[2] += (CD_ca_atoms[ci].position[2])/np.float(len(CD_ca_atoms))
        
        for ei in range(len(EF_ca_atoms)):
            EF_com[0] += (EF_ca_atoms[ei].position[0])/np.float(len(EF_ca_atoms))
            EF_com[1] += (EF_ca_atoms[ei].position[1])/np.float(len(EF_ca_atoms))
            EF_com[2] += (EF_ca_atoms[ei].position[2])/np.float(len(EF_ca_atoms))
            
        for gi in range(len(GH_ca_atoms)):
            GH_com[0] += (GH_ca_atoms[gi].position[0])/np.float(len(GH_ca_atoms))
            GH_com[1] += (GH_ca_atoms[gi].position[1])/np.float(len(GH_ca_atoms))
            GH_com[2] += (GH_ca_atoms[gi].position[2])/np.float(len(GH_ca_atoms))
        
        for ii in range(len(IJ_ca_atoms)):
            IJ_com[0] += (IJ_ca_atoms[ii].position[0])/np.float(len(IJ_ca_atoms))
            IJ_com[1] += (IJ_ca_atoms[ii].position[1])/np.float(len(IJ_ca_atoms))
            IJ_com[2] += (IJ_ca_atoms[ii].position[2])/np.float(len(IJ_ca_atoms))
        #alpha_array.append(angle_between_3(AB_com[0], AB_com[1], AB_com[2], tot_com[0], tot_com[1], tot_com[2], EF_com[0], EF_com[1], EF_com[2]))
        #gamma_array.append(angle_between_3(CD_com[0], CD_com[1], CD_com[2], tot_com[0], tot_com[1], tot_com[2], GH_com[0], GH_com[1], GH_com[2]))
        #beta1_array.append(angle_between_3(AB_com[0], AB_com[1], AB_com[2], tot_com[0], tot_com[1], tot_com[2], CD_com[0], CD_com[1], CD_com[2]))
        #beta2_array.append(angle_between_3(EF_com[0], EF_com[1], EF_com[2], tot_com[0], tot_com[1], tot_com[2], GH_com[0], GH_com[1], GH_com[2]))

        AB_to_CD.append(distance_between(AB_com[0], AB_com[1], AB_com[2],CD_com[0], CD_com[1], CD_com[2]))
        AB_to_EF.append(distance_between(AB_com[0], AB_com[1], AB_com[2],EF_com[0], EF_com[1], EF_com[2]))
        EF_to_GH.append(distance_between(EF_com[0], EF_com[1], EF_com[2],GH_com[0], GH_com[1], GH_com[2]))
        CD_to_GH.append(distance_between(CD_com[0], CD_com[1], CD_com[2],GH_com[0], GH_com[1], GH_com[2]))
        IJ_to_AB.append(distance_between(IJ_com[0], IJ_com[1], IJ_com[2],AB_com[0], AB_com[1], AB_com[2]))
        IJ_to_GH.append(distance_between(IJ_com[0], IJ_com[1], IJ_com[2],GH_com[0], GH_com[1], GH_com[2]))

        for r in residue_array:
            r_string = "resid " + str(r)
            r_atoms = u.select_atoms(r_string)
            r_length = len(r_atoms)
            for n in range(r_length):
                r_atoms[n].position[0] -= multiplier_list[m]*pc_1[3*(r-1)]
                r_atoms[n].position[1] -= multiplier_list[m]*pc_1[3*(r-1)+1]
                r_atoms[n].position[2] -= multiplier_list[m]*pc_1[3*(r-1)+2]

w.close_trajectory()

"""begin
plt.figure()
plt.xlabel(r"PC Multiplier")
plt.ylabel(r"Distance between Dimers ($\AA$)")
plt.title("H3 NUC")
plt.xlim((-1.0*multiplier,1.0*multiplier))
plt.ylim((31.0,38.0))
plt.plot(multiplier_array_a,AB_to_EF, color="blue", linestyle="-", linewidth=1.5, label = "ABEF homotetramer")
plt.plot(multiplier_array_a,CD_to_GH, color="green", linestyle="-", linewidth=1.5, label = "CDGH homotetramer")
plt.plot(multiplier_array_a,AB_to_CD, color="orange", linestyle="-", linewidth=1.5, label = "ABCD heterotetramer")
plt.plot(multiplier_array_a,EF_to_GH, color="orangered", linestyle="-", linewidth=1.5, label = "EFGH heterotetramer")
plt.axvline(x=0.0, color="gray", linestyle="--", linewidth=0.5)
end"""
#legend0 = plt.legend(loc='upper right')
#frame0  = legend0.get_frame()
#frame0.set_facecolor('1.00')
#for label in legend0.get_texts():
#    label.set_fontsize(8)
#plt.show()


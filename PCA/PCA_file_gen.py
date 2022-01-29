from math import *
from numpy import *
from matplotlib import *
from pylab import *
import MDAnalysis
import numpy
import numpy.linalg
import matplotlib.pyplot as plt
import io


u = MDAnalysis.Universe(GRO,TRJ)
cut = 0 # times 25(frames skipped)*2(ps/frame) = 600000 ps cut                                                     
start = cut + 1
timestep = 1
n_frames = len(u.trajectory)-1						
print "n_frames: ", n_frames

# Currenlt ca_atoms is set up for you to enter the amount, this can be automated. 
ca_atoms = u.selectAtoms('bynum 213:2787 and name CA') + u.selectAtoms('bynum 2947:8501 and name CA') + u.selectAtoms('bynum 8642:10058 and name CA') + u.selectAtoms('bynum 10240:11552 and name CA') + u.selectAtoms('bynum 12337:15389 and name P') + u.selectAtoms('bynum 16979:20031 and name P')
print "len(ca_atoms): ", len(ca_atoms)
print ca_atoms

# calculate com_total
tot_com = numpy.zeros((n_frames,3))					
for ts in u.trajectory[start::timestep]:                                                                                                                                                                   
    n = ts.frame - 2                                                                                                                                                                                  
    print n
    for ti in range(len(ca_atoms)):
        tot_com[n][0] += ca_atoms[ti].pos[0]                                                                                                                                                          
        tot_com[n][1] += ca_atoms[ti].pos[1]                                                                                                                                                            
        tot_com[n][2] += ca_atoms[ti].pos[2]                                                                                                                                                            
    tot_com[n][0] /= numpy.float(len(ca_atoms))                                                                                                                                                       
    tot_com[n][1] /= numpy.float(len(ca_atoms))                                                                                                                                                        
    tot_com[n][2] /= numpy.float(len(ca_atoms))   


# calculate PCA input file, CA positions with respect to COM
f_ca_atoms = open('p_ca_t.output','w')                      
for ts in u.trajectory[start::timestep]:
    t = ts.frame - 2
    print "t:", t
    for i in range(len(ca_atoms)):
	print "i:", i
        ca_x = ca_atoms[i].pos[0] - tot_com[t][0]
        ca_y = ca_atoms[i].pos[1] - tot_com[t][1]
        ca_z = ca_atoms[i].pos[2] - tot_com[t][2]
        position_string = str(ca_x) + ' ' + str(ca_y) + ' ' + str(ca_z) + ' '
        f_ca_atoms.write(position_string)
    f_ca_atoms.write('\n')
f_ca_atoms.close()

"""beginning

for i in range(len(ca_atoms)):
    with io.open("position_ca_" + str(i) + ".dat", 'w', encoding='utf-8') as f:
       f.write(str(func(i))



A_ca_atoms = ca_atoms[0:89]
B_ca_atoms = ca_atoms[89:167]
C_ca_atoms = ca_atoms[167:264]
D_ca_atoms = ca_atoms[264:354]
E_ca_atoms = ca_atoms[354:442]
F_ca_atoms = ca_atoms[442:521]
G_ca_atoms = ca_atoms[521:621]
H_ca_atoms = ca_atoms[621:711]
#print A_ca_atoms
#print A_ca_atoms[0]
#print A_ca_atoms[88]
#print B_ca_atoms[0]
#print len(A_ca_atoms), len(A_ca_atoms)+len(B_ca_atoms), len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms), len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms), len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms)+len(E_ca_atoms), len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms)+len(E_ca_atoms)+len(F_ca_atoms),len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms)+len(E_ca_atoms)+len(F_ca_atoms)+len(G_ca_atoms),len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms)+len(E_ca_atoms)+len(F_ca_atoms)+len(G_ca_atoms)+len(H_ca_atoms)

total_length = len(A_ca_atoms)+len(B_ca_atoms)+len(C_ca_atoms)+len(D_ca_atoms)+len(E_ca_atoms)+len(F_ca_atoms)+len(G_ca_atoms)+len(H_ca_atoms)

print "check: ",(total_length==711)

tot_com = numpy.zeros((n_frames,3))
A_com = numpy.zeros((n_frames,3))
B_com = numpy.zeros((n_frames,3))
C_com = numpy.zeros((n_frames,3))
D_com = numpy.zeros((n_frames,3))
E_com = numpy.zeros((n_frames,3))
F_com = numpy.zeros((n_frames,3))
G_com = numpy.zeros((n_frames,3))
H_com = numpy.zeros((n_frames,3))


# calculate center-of-mass matrices
for ts in u.trajectory[0::timestep]:
    n = ts.frame - 1
    print n
    
    for ti in range(len(ca_atoms)):
        tot_com[n][0] += ca_atoms[ti].pos[0]
        tot_com[n][1] += ca_atoms[ti].pos[1]
        tot_com[n][2] += ca_atoms[ti].pos[2]
    tot_com[n][0] /= numpy.float(len(ca_atoms))
    tot_com[n][1] /= numpy.float(len(ca_atoms))
    tot_com[n][2] /= numpy.float(len(ca_atoms))

    for ai in range(len(A_ca_atoms)):
        A_com[n][0] += A_ca_atoms[ai].pos[0]
        A_com[n][1] += A_ca_atoms[ai].pos[1]
        A_com[n][2] += A_ca_atoms[ai].pos[2]
    A_com[n][0] /= numpy.float(len(A_ca_atoms))
    A_com[n][1] /= numpy.float(len(A_ca_atoms))
    A_com[n][2] /= numpy.float(len(A_ca_atoms))
    
    for bi in range(len(B_ca_atoms)):
        B_com[n][0] += B_ca_atoms[bi].pos[0]
        B_com[n][1] += B_ca_atoms[bi].pos[1]
        B_com[n][2] += B_ca_atoms[bi].pos[2]
    B_com[n][0] /= numpy.float(len(B_ca_atoms))
    B_com[n][1] /= numpy.float(len(B_ca_atoms))
    B_com[n][2] /= numpy.float(len(B_ca_atoms))

    for ci in range(len(C_ca_atoms)):
        C_com[n][0] += C_ca_atoms[ci].pos[0]
        C_com[n][1] += C_ca_atoms[ci].pos[1]
        C_com[n][2] += C_ca_atoms[ci].pos[2]
    C_com[n][0] /= numpy.float(len(C_ca_atoms))
    C_com[n][1] /= numpy.float(len(C_ca_atoms))
    C_com[n][2] /= numpy.float(len(C_ca_atoms))

    for di in range(len(D_ca_atoms)):
        D_com[n][0] += D_ca_atoms[di].pos[0]
        D_com[n][1] += D_ca_atoms[di].pos[1]
        D_com[n][2] += D_ca_atoms[di].pos[2]
    D_com[n][0] /= numpy.float(len(D_ca_atoms))
    D_com[n][1] /= numpy.float(len(D_ca_atoms))
    D_com[n][2] /= numpy.float(len(D_ca_atoms))

    for ei in range(len(E_ca_atoms)):
        E_com[n][0] += E_ca_atoms[ei].pos[0]
        E_com[n][1] += E_ca_atoms[ei].pos[1]
        E_com[n][2] += E_ca_atoms[ei].pos[2]
    E_com[n][0] /= numpy.float(len(E_ca_atoms))
    E_com[n][1] /= numpy.float(len(E_ca_atoms))
    E_com[n][2] /= numpy.float(len(E_ca_atoms))

    for fi in range(len(F_ca_atoms)):
        F_com[n][0] += F_ca_atoms[fi].pos[0]
        F_com[n][1] += F_ca_atoms[fi].pos[1]
        F_com[n][2] += F_ca_atoms[fi].pos[2]
    F_com[n][0] /= numpy.float(len(F_ca_atoms))
    F_com[n][1] /= numpy.float(len(F_ca_atoms))
    F_com[n][2] /= numpy.float(len(F_ca_atoms))

    for gi in range(len(G_ca_atoms)):
        G_com[n][0] += G_ca_atoms[gi].pos[0]
        G_com[n][1] += G_ca_atoms[gi].pos[1]
        G_com[n][2] += G_ca_atoms[gi].pos[2]
    G_com[n][0] /= numpy.float(len(G_ca_atoms))
    G_com[n][1] /= numpy.float(len(G_ca_atoms))
    G_com[n][2] /= numpy.float(len(G_ca_atoms))

    for hi in range(len(H_ca_atoms)):
        H_com[n][0] += H_ca_atoms[hi].pos[0]
        H_com[n][1] += H_ca_atoms[hi].pos[1]
        H_com[n][2] += H_ca_atoms[hi].pos[2]
    H_com[n][0] /= numpy.float(len(H_ca_atoms))
    H_com[n][1] /= numpy.float(len(H_ca_atoms))
    H_com[n][2] /= numpy.float(len(H_ca_atoms))

# calculate d(COM[i],COM[j]), distances between centers-of-mass

f_00 = open('d_COM_00.output','w')
f_01 = open('d_COM_01.output','w')
f_02 = open('d_COM_02.output','w')
f_03 = open('d_COM_03.output','w')
f_04 = open('d_COM_04.output','w')
f_05 = open('d_COM_05.output','w')
f_06 = open('d_COM_06.output','w')
f_07 = open('d_COM_07.output','w')
f_08 = open('d_COM_08.output','w')
f_09 = open('d_COM_09.output','w')
f_10 = open('d_COM_10.output','w')
f_11 = open('d_COM_11.output','w')
f_12 = open('d_COM_12.output','w')
f_13 = open('d_COM_13.output','w')
f_14 = open('d_COM_14.output','w')
f_15 = open('d_COM_15.output','w')
f_16 = open('d_COM_16.output','w')
f_17 = open('d_COM_17.output','w')
f_18 = open('d_COM_18.output','w')
f_19 = open('d_COM_19.output','w')
f_20 = open('d_COM_20.output','w')
f_21 = open('d_COM_21.output','w')
f_22 = open('d_COM_22.output','w')
f_23 = open('d_COM_23.output','w')

d_COM = numpy.zeros((n_frames,24))

for t in range(n_frames):
    print t
    
    
    Ax_com = A_com[t][0]-tot_com[t][0]
    Ay_com = A_com[t][1]-tot_com[t][1]
    Az_com = A_com[t][2]-tot_com[t][2]
    f_00.write(str(t) + ' ' + str(Ax_com) + '\n')
    f_01.write(str(t) + ' ' + str(Ay_com) + '\n')
    f_02.write(str(t) + ' ' + str(Az_com) + '\n')


    Bx_com = B_com[t][0]-tot_com[t][0]
    By_com = B_com[t][1]-tot_com[t][1]
    Bz_com = B_com[t][2]-tot_com[t][2]
    f_03.write(str(t) + ' ' + str(Bx_com) + '\n')
    f_04.write(str(t) + ' ' + str(By_com) + '\n')
    f_05.write(str(t) + ' ' + str(Bz_com) + '\n')
    

    Cx_com = C_com[t][0]-tot_com[t][0]
    Cy_com = C_com[t][1]-tot_com[t][1]
    Cz_com = C_com[t][2]-tot_com[t][2]
    f_06.write(str(t) + ' ' + str(Cx_com) + '\n')
    f_07.write(str(t) + ' ' + str(Cy_com) + '\n')
    f_08.write(str(t) + ' ' + str(Cz_com) + '\n')


    Dx_com = D_com[t][0]-tot_com[t][0]
    Dy_com = D_com[t][1]-tot_com[t][1]
    Dz_com = D_com[t][2]-tot_com[t][2]
    f_09.write(str(t) + ' ' + str(Dx_com) + '\n')
    f_10.write(str(t) + ' ' + str(Dy_com) + '\n')
    f_11.write(str(t) + ' ' + str(Dz_com) + '\n')


    Ex_com = E_com[t][0]-tot_com[t][0]
    Ey_com = E_com[t][1]-tot_com[t][1]
    Ez_com = E_com[t][2]-tot_com[t][2]
    f_12.write(str(t) + ' ' + str(Ex_com) + '\n')
    f_13.write(str(t) + ' ' + str(Ey_com) + '\n')
    f_14.write(str(t) + ' ' + str(Ez_com) + '\n')


    Fx_com = F_com[t][0]-tot_com[t][0]
    Fy_com = F_com[t][1]-tot_com[t][1]
    Fz_com = F_com[t][2]-tot_com[t][2]
    f_15.write(str(t) + ' ' + str(Fx_com) + '\n')
    f_16.write(str(t) + ' ' + str(Fy_com) + '\n')
    f_17.write(str(t) + ' ' + str(Fz_com) + '\n')

    Gx_com = G_com[t][0]-tot_com[t][0]
    Gy_com = G_com[t][1]-tot_com[t][1]
    Gz_com = G_com[t][2]-tot_com[t][2]
    f_18.write(str(t) + ' ' + str(Gx_com) + '\n')
    f_19.write(str(t) + ' ' + str(Gy_com) + '\n')
    f_20.write(str(t) + ' ' + str(Gz_com) + '\n')

    Hx_com = H_com[t][0]-tot_com[t][0]
    Hy_com = H_com[t][1]-tot_com[t][1]
    Hz_com = H_com[t][2]-tot_com[t][2]
    f_21.write(str(t) + ' ' + str(Hx_com) + '\n')
    f_22.write(str(t) + ' ' + str(Hy_com) + '\n')
    f_23.write(str(t) + ' ' + str(Hz_com) + '\n')
    
# check individual files after this script runs
end"""

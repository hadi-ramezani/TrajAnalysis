#! /usr/bin/env python

from MDAnalysis import *
from math import *
import numpy as np
import time
import sys
import argparse
import os
import matplotlib.pyplot as plt
from subprocess import *
from pylab import *

parser = argparse.ArgumentParser(description="Read trajectory and coordinate file")

parser.add_argument("-t",
                    action="store", nargs='?', default="",
                    required=False, dest="traj",
                    help="specifies a dcd file created using NAMD package")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-d",
                    action="store", nargs='?',
                    required=True, dest="dz", type=float,
                    help="Grid thickness")
parser.add_argument("-ti",
                    action="store", nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/wat_dielectric.dat",
                    required=False, dest="output",
                    help="output filename")
parser.add_argument("-i",
                    action="store", nargs='?', default="",
                    required=False, dest="input",
                    help="input filename")
parser.add_argument("-b",
                    action="store", nargs='?', default=int('0'),
                    required=False, dest="start_time", type=int,
                    help="First time to start reading trajectory")
parser.add_argument("-e",
                    action="store", nargs='?', default=float('inf'),
                    required=False, dest="end_time", type=int,
                    help="Last time to stop readining trajectory")
parser.add_argument("-f",
                    action="store", nargs='?',
                    required=False, dest="firstFileName", type=int,
                    help="The name of the first trajectory file before dot")
parser.add_argument("-l",
                    action="store", nargs='?',
                    required=False, dest="lastFileName", type=int,
                    help="The name of the last trajectory file before dot")
parser.add_argument("-s",
                    action="store", nargs='?', default="name OH2",
                    required=False, dest="atom_selection",
                    help="Atom selection argument")

args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
dz = args.dz
output_filename = args.output
input_filename = args.input
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection = args.atom_selection

time1 = time.time()
output_file = open(output_filename, "w")

# Define Constants
#k = 1.240865258e-23
coef = 332.0636
k= 1.987204395e-3
T = 298.15
kT = (k*T)/coef
kTinv = pow(kT, -1)

def vector21(v1, v2):
    return [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]]
def position(v1, v2):
    return (v1 + v2) / 2.0

# Read the trj and psf file
if traj_filename == "":
    list = [str(firstFile) + ".dcd"]
    for i in xrange(firstFile + 1, lastFile + 1):
        list.append(str(i) + ".dcd")
    print "List of Trajectory files:" + str(list)
    u = Universe(psf_filename, list)
else:
    u = Universe(psf_filename, traj_filename)

# Construct the B vector for the selected phase
ind = u.selectAtoms(atom_selection).residues[0].atoms.indices()
q = u.selectAtoms(atom_selection).residues[0].charges()
bond = u.selectAtoms(atom_selection).residues[0].bonds.dump_contents()

C = np.empty((len(q), len(bond)), dtype='float')

for i in xrange(0, len(ind)):
    atomi = ind[i]
    for l in xrange(0, len(bond)):
        atomj = bond[l][0]
        atomk = bond[l][1]
        if atomi == atomj:
            C[i][l] = 1
        elif atomi == atomk:
            C[i][l] = -1
        else:
            C[i][l] = 0

C_inv = pinv(C)
q = q.reshape((len(q), 1))
B = dot(C_inv, q)
len_B = len(B)

# Determine the number of molecules
n_mol = u.selectAtoms(atom_selection)

# Get the atoms indices of bonds for molecules
bond_all = u.selectAtoms(atom_selection).residues.bonds.dump_contents()
bond_all = np.array(bond_all)

# Initialize matrices to store coordinates and dipole-moments
coord = np.empty((len(bond_all), 6))
mu = np.empty((len(bond_all), 6))


# Obtain initial information form psf and trj files
natoms = len(u.selectAtoms("all"))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)

# Get number of frames
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames)
if end_frame > num_frames:
    end_frame = num_frames
print "reading frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays (the ones which are independent of nslabs)
pvector = np.zeros((natoms), dtype='float')

P_f_x = np.empty((0), dtype='float')
P_f_y = np.empty((0), dtype='float')
P_f_xy = np.empty((0), dtype='float')
P_f_z = np.empty((0), dtype='float')
M_f_x = np.empty((0), dtype='float')
M_f_y = np.empty((0), dtype='float')
M_f_z = np.empty((0), dtype='float')
M2_f_x = np.empty((0), dtype='float')
M2_f_y = np.empty((0), dtype='float')
M2_f_z = np.empty((0), dtype='float')
MP_f_x = np.empty((0), dtype='float')
MP_f_y = np.empty((0), dtype='float')
MP_f_z = np.empty((0), dtype='float')

nframes = 0

# Get the box Size Information
box = u.dimensions[0:3]

# Loop over selected frames
for curr_frame in xrange(start_frame, end_frame) :
    if curr_frame != 0 :
        trj = u.trajectory[curr_frame]
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame >= start_frame and curr_frame <= end_frame :
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
        # Get coordinates of all atoms (to be used to retrieve the coordinates of molecules)
        coord_all = u.selectAtoms("all").residues.coordinates()
        coord_all = np.array(coord_all)

        # Get the coordinates and calculate dipole-moments
        for i in xrange(0, len(n_mol)):
            for j in xrange(0, len_B):
                index = i * len_B + j
                coord[index][0:3] = coord_all[bond_all[index][0]]
                coord[index][3:6] = coord_all[bond_all[index][1]]
                mu[index][0:3] = position(coord[index][0:3], coord[index][3:6])
                mu[index][3:6] = B[j] * vector21(coord[index][0:3], coord[index][3:6])


        M_x = np.sum(mu[:, 3])
        M_y = np.sum(mu[:, 4])
        M_z = np.sum(mu[:, 5])

        M2_x = np.power(M_x, 2)
        M2_y = np.power(M_y, 2)
        M2_z = np.power(M_z, 2)

        M_f_x = np.append(M_f_x, M_x)
        M_f_y = np.append(M_f_y, M_y)
        M_f_z = np.append(M_f_z, M_z)

        M2_f_x = np.append(M2_f_x, M2_x)
        M2_f_y = np.append(M2_f_y, M2_y)
        M2_f_z = np.append(M2_f_z, M2_z)

        # Define slabs
        if curr_frame == start_frame:
            for i in xrange(0, len(coord_all)):
                loc = coord_all[i]
                pvector[i] = loc[2]
            pvector = np.sort(pvector)
            lowZ = int(pvector[0])
            highZ = int(pvector[-1])
            bins = np.arange(lowZ - 3*dz, highZ + 3*dz, dz)
            nslab = len(bins) - 1
            # Initialize arrays
            mu_x = np.empty((nslab), dtype='float')
            mu_y = np.empty((nslab), dtype='float')
            mu_z = np.empty((nslab), dtype='float')
            P_x = np.empty((nslab), dtype='float')
            P_y = np.empty((nslab), dtype='float')
            P_z = np.empty((nslab), dtype='float')


        mu_x.fill(0)
        mu_y.fill(0)
        mu_z.fill(0)
        P_x.fill(0)
        P_y.fill(0)
        P_z.fill(0)

        for i in xrange(0, len(bond_all)):
            p = mu[i][2]
            x = int((p-bins[0])/dz)

            if x < nslab and x >= 0:
                mu_x[x] += mu[i][3]
                mu_y[x] += mu[i][4]
                mu_z[x] += mu[i][5]


        for i in xrange(nslab):
            P_x[i] = mu_x[i] / box[0] / box[1] / dz
            P_y[i] = mu_y[i] / box[0] / box[1] / dz
            P_z[i] = mu_z[i] / box[0] / box[1] / dz

        P_f_z = np.append(P_f_z, P_z, 0)
        P_f_x = np.append(P_f_x, P_x, 0)
        P_f_y = np.append(P_f_y, P_y, 0)

        MP_f_x = np.append(MP_f_x, P_x * M_x, 0)
        MP_f_y = np.append(MP_f_y, P_y * M_y, 0)
        MP_f_z = np.append(MP_f_z, P_z * M_z, 0)

P_f_x = np.reshape(P_f_x, (len(P_f_x) / nslab, nslab))
P_f_x_ave = np.average(P_f_x, 0)
P_f_x_dev = np.std(P_f_x, 0)

P_f_y = np.reshape(P_f_y, (len(P_f_y) / nslab, nslab))
P_f_y_ave = np.average(P_f_y, 0)
P_f_y_dev = np.std(P_f_y, 0)

P_f_z = np.reshape(P_f_z, (len(P_f_z) / nslab, nslab))
P_f_z_ave = np.average(P_f_z, 0)
P_f_z_dev = np.std(P_f_z, 0)

MP_f_x = np.reshape(MP_f_x, (len(MP_f_x) / nslab, nslab))
MP_f_x_ave = np.average(MP_f_x, 0)
MP_f_x_dev = np.std(MP_f_x, 0)

MP_f_y = np.reshape(MP_f_y, (len(MP_f_y) / nslab, nslab))
MP_f_y_ave = np.average(MP_f_y, 0)
MP_f_y_dev = np.std(MP_f_y, 0)

MP_f_z = np.reshape(MP_f_z, (len(MP_f_z) / nslab, nslab))
MP_f_z_ave = np.average(MP_f_z, 0)
MP_f_z_dev = np.std(MP_f_z, 0)

M_f_x_ave = np.average(M_f_x)
M_f_y_ave = np.average(M_f_y)
M_f_z_ave = np.average(M_f_z)

M2_f_x_ave = np.average(M2_f_x)
M2_f_y_ave = np.average(M2_f_y)
M2_f_z_ave = np.average(M2_f_z)

# Define h, H, and e (dielectric constant)
h_xy = np.empty((nslab),dtype='float')
h_z  = np.empty((nslab),dtype='float')

X_xy = np.empty((nslab),dtype='float')
X_z  = np.empty((nslab),dtype='float')

e_xy = np.empty((nslab),dtype='float')
e_z  = np.empty((nslab),dtype='float')

h_xy = 0.5 * kTinv * ((MP_f_x_ave - P_f_x_ave * M_f_x_ave) + (MP_f_y_ave - P_f_y_ave * M_f_y_ave))
h_z = kTinv * (MP_f_z_ave - P_f_z_ave * M_f_z_ave)
H_xy = 0.5 * kTinv * ((M2_f_x_ave - M_f_x_ave * M_f_x_ave)+(M2_f_y_ave - M_f_y_ave * M_f_y_ave))
H_z = kTinv * (M2_f_z_ave - M_f_z_ave * M_f_z_ave)

V = box[0] * box [1] * box[2]

X_xy = h_xy
X_z = h_z * np.power(1+ 4 * pi *(H_z/V - h_z), -1)

# Get the values of dielectric constant
e_xy = 4 * pi * X_xy + 1
e_z = 4 * pi * X_z + 1

# Write the data into a file
for x in range(nslab):
    p = (bins[x] + bins[x + 1]) / 2
    output_file.write('%5.3f %10.7f %10.7f\n' % (
    p, e_z[x], e_xy[x]))

time2 = time.time()
S_time = time2 - time1
print "Simulation time = " + str(S_time)
output_file.close()

# Plot the components of dielectric constant
def read_data(output_filename):
    outFile = open(str(output_filename), 'r')
    binCenter = []
    e_z = []
    e_xy = []

    for line in outFile:
        line = line.strip()
        if line != "":
            columns = line.split()
            if '#' in columns[0]:
                check = 'yes'
            else:
                firstCol = float(columns[0])
                binCenter.append(firstCol)
                secCol = float(columns[1])
                e_z.append(secCol)
                thirCol = float(columns[2])
                e_xy.append(thirCol)

    e_z = np.array(e_z)
    e_xy = np.array(e_xy)
    return binCenter, e_z, e_xy

def plot_data(binCenter, e_xy, color, label):
    plt.rcParams['axes.linewidth'] = 1.5
    # plt.xticks(rotation=25)
    plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
    plt.ylabel(r'$\mathrm{\epsilon_{=}(z)}$', fontsize=23, labelpad=20)
    plt.tick_params(which='both', width=2)
    plt.xlim((binCenter[0],binCenter[len(binCenter)-1]))
    #plt.ylim((-40, 140))
    plt.yticks(fontsize=17)
    plt.xticks(fontsize=17)
    plt.gcf().set_size_inches(10, 5)
    plt.plot(binCenter, e_xy, color=color, linewidth=2, label=label)
    plt.legend(loc='best', frameon=True, fancybox=True)

binCenter, e_z, e_xy = read_data(output_filename)
plot_data(binCenter, e_xy, 'red', "Water")

if input_filename != "":
    binCenter, e_z, e_xy = read_data(input_filename)
    plot_data(binCenter, e_xy, 'blue', "5CB")

plt.savefig('../Analysis/plot_dielectric.png', dpi=600, facecolor='w', edgecolor='w',
            orientation='portrait', papertype='letter', format='png',
            transparent=True, bbox_inches='tight', pad_inches=0.1,
            frameon=None)
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
                    action="store", nargs='?', default="../Analysis/Polarization.dat",
                    required=False, dest="output",
                    help="output filename")
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

args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
dz = args.dz
output_filename = args.output
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection_N = "name NY1"
atom_selection_C = "name CA12"

time1 = time.time()
output_file = open(output_filename, "w")


def distance(v1, v2):
    return sqrt(
        (v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) + (v1[2] - v2[2]) * (v1[2] - v2[2]))
def distance_d(v1, v2):
    return abs(v1 - v2)
def vector12(v1, v2):
    return [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]
def vector21(v1, v2):
    return [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]]
def position(v1, v2):
    return (v1 + v2) / 2.0
def Sz(v):
    return (3. * v[2] * v[2] - 1.) / 2.0


# Read the trj and psf file
if traj_filename == "":
    list = [str(firstFile) + ".dcd"]
    for i in xrange(firstFile + 1, lastFile + 1):
        list.append(str(i) + ".dcd")
    print "List of Trajectory files:" + str(list)
    u = Universe(psf_filename, list)
else:
    u = Universe(psf_filename, traj_filename)

# Construct the B vector for water phase
ind_wat = u.selectAtoms("name OH2").residues[0].atoms.indices()
q_wat = u.selectAtoms("name OH2").residues[0].charges()
bond_wat = u.selectAtoms("name OH2").residues[0].bonds.dump_contents()

C_wat = np.empty((len(q_wat), len(bond_wat)), dtype='float')

for i in xrange(0, len(ind_wat)):
    atomi = ind_wat[i]
    for l in xrange(0, len(bond_wat)):
        atomj = bond_wat[l][0]
        atomk = bond_wat[l][1]
        if atomi == atomj:
            C_wat[i][l] = 1
        elif atomi == atomk:
            C_wat[i][l] = -1
        else:
            C_wat[i][l] = 0

C_wat_inv = pinv(C_wat)
q_wat = q_wat.reshape((len(q_wat), 1))
B_wat = dot(C_wat_inv, q_wat)
len_B_wat = len(B_wat)

# Construct the B vector for LC phase
ind_lc = u.selectAtoms("name NY1").residues[0].atoms.indices()
q_lc = u.selectAtoms("name NY1").residues[0].charges()
bond_lc = u.selectAtoms("name NY1").residues[0].bonds.dump_contents()

C_lc = np.empty((len(q_lc), len(bond_lc)), dtype='float')

for i in xrange(0, len(ind_lc)):
    atomi = ind_lc[i]
    for l in xrange(0, len(bond_lc)):
        atomj = bond_lc[l][0]
        atomk = bond_lc[l][1]
        if atomi == atomj:
            C_lc[i][l] = 1
        elif atomi == atomk:
            C_lc[i][l] = -1
        else:
            C_lc[i][l] = 0

C_lc_inv = pinv(C_lc)
q_lc = q_lc.reshape((len(q_lc), 1))
B_lc = dot(C_lc_inv, q_lc)
len_B_lc = len(B_lc)

# Determine the number of water molecules and the number of liquid crystal molecules
nwat = u.selectAtoms("name OH2")
nlc = u.selectAtoms("name NY1")

# Get the atoms indices of bonds for water and LC molecules
bond_wat_all = u.selectAtoms("name OH2").residues.bonds.dump_contents()
bond_wat_all = np.array(bond_wat_all)
bond_lc_all = u.selectAtoms("name NY1").residues.bonds.dump_contents()
bond_lc_all = np.array(bond_lc_all)

# Initialize matrices to store coordinates and dipole-moments
wat_coord = np.empty((len(bond_wat_all), 6))
mu_wat = np.empty((len(bond_wat_all), 6))
lc_coord = np.empty((len(bond_lc_all), 6))
mu_lc = np.empty((len(bond_lc_all), 6))


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

# Deine arrays for water and LC (the ones which are independent of nslabs)
pvector = np.zeros((natoms), dtype='float')

wat_P_Z_f_x = np.empty((0), dtype='float')
wat_P_Z_f_y = np.empty((0), dtype='float')
wat_P_Z_f_xy = np.empty((0), dtype='float')
wat_P_Z_f_z = np.empty((0), dtype='float')

lc_P_Z_f_x = np.empty((0), dtype='float')
lc_P_Z_f_y = np.empty((0), dtype='float')
lc_P_Z_f_xy = np.empty((0), dtype='float')
lc_P_Z_f_z = np.empty((0), dtype='float')

# Loop over selected frames
nframes = 0

# Get the box Size Information
box = u.dimensions[0:3]

for curr_frame in xrange(0, num_frames) :
    if curr_frame != 0 :
        trj = u.trajectory[curr_frame]
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame > end_frame :
        break
    else :
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
        # Get coordinates of all atoms (to be used to retrieve the coordinates of water and lc molecules)
        coord_all = u.selectAtoms("all").residues.coordinates()
        coord_all = np.array(coord_all)

        # Get the coordinates and calculate dipole-moments of water phase
        for i in xrange(0, len(nwat)):
            for j in xrange(0, len_B_wat):
                index = i * len_B_wat + j
                wat_coord[index][0:3] = coord_all[bond_wat_all[index][0]]
                wat_coord[index][3:6] = coord_all[bond_wat_all[index][1]]
                mu_wat[index][0:3] = position(wat_coord[index][0:3],
                                                           wat_coord[index][3:6])
                mu_wat[index][3:6] = B_wat[j] * vector21(wat_coord[index][0:3],
                                                                      wat_coord[index][3:6])
        # Get the coordinates and calculate dipole-moments of lc phase
        for i in xrange(0, len(nlc)):
            for j in xrange(0, len_B_lc):
                index = i * len_B_lc + j
                lc_coord[index][0:3] = coord_all[bond_lc_all[index][0]]
                lc_coord[index][3:6] = coord_all[bond_lc_all[index][1]]
                mu_lc[index][0:3] = position(lc_coord[index][0:3],
                                                         lc_coord[index][3:6])
                mu_lc[index][3:6] = B_lc[j] * vector21(lc_coord[index][0:3],
                                                                   lc_coord[index][3:6])

        # Define slabs
        if curr_frame == start_frame:
            for i in xrange(0, len(coord_all)):
                loc = coord_all[i]
                pvector[i] = loc[2]
            pvector = np.sort(pvector)
            highZ = int(pvector[natoms - 1]) + 5
            lowZ = int(pvector[0])
            nslab = int((highZ - lowZ) / dz)
            bins = np.zeros((nslab + 1), dtype='float')
            # Initialize arrayes
            wat_M_Z_z = np.empty((nslab), dtype='float')
            lc_M_Z_z = np.empty((nslab), dtype='float')

            wat_P_Z_z = np.empty((nslab), dtype='float')
            lc_P_Z_z = np.empty((nslab), dtype='float')

            bins[nslab] = highZ
            bins[0] = lowZ
            slab_ind = nslab -1
            for i in xrange(1, nslab):
                bins[i] = bins[i - 1] + dz

        wat_M_Z_z.fill(0)
        lc_M_Z_z.fill(0)

        wat_P_Z_z.fill(0)
        lc_P_Z_z.fill(0)

        for i in xrange(0, len(bond_wat_all)):
            wat_p = mu_wat[i][2]
            x = min(int((wat_p-bins[0])/dz),slab_ind)
            wat_M_Z_z[x] += mu_wat[i][5]

        for i in xrange(0, len(bond_lc_all)):
            lc_p = mu_lc[i][2]
            x = min(int((lc_p-bins[0])/dz),slab_ind)
            lc_M_Z_z[x] += mu_lc[i][5]

        for i in xrange(nslab):
            wat_P_Z_z[i] = wat_M_Z_z[i] / box[0] / box[1] / dz
            lc_P_Z_z[i] = lc_M_Z_z[i] / box[0] / box[1] / dz

        wat_P_Z_f_z = np.append(wat_P_Z_f_z, wat_P_Z_z, 0)
        lc_P_Z_f_z = np.append(lc_P_Z_f_z, lc_P_Z_z, 0)

wat_P_Z_f_z = np.reshape(wat_P_Z_f_z, (len(wat_P_Z_f_z) / nslab, nslab))
lc_P_Z_f_z = np.reshape(lc_P_Z_f_z, (len(lc_P_Z_f_z) / nslab, nslab))

wat_P_Z_f_z_ave = np.average(wat_P_Z_f_z, 0)
lc_P_Z_f_z_ave = np.average(lc_P_Z_f_z, 0)

wat_P_Z_f_z_dev = np.std(wat_P_Z_f_z, 0)
lc_P_Z_f_z_dev = np.std(lc_P_Z_f_z, 0)

for x in range(nslab):
    p = (bins[x] + bins[x + 1]) / 2
    output_file.write('%5.3f %10.7f %10.7f %10.7f %10.7f\n' % (
    p, wat_P_Z_f_z_ave[x], wat_P_Z_f_z_dev[x], lc_P_Z_f_z_ave[x], lc_P_Z_f_z_dev[x]))

time2 = time.time()
S_time = time2 - time1
print "Simulation time = " + str(S_time)
output_file.close()

# Plot the components of polarization
outFile = open(str(output_filename), 'r')
binCenter = []
wat_Pz = []
wat_Pz_dev = []

lc_Pz = []
lc_Pz_dev = []

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
            wat_Pz.append(secCol)
            thirCol = float(columns[2])
            wat_Pz_dev.append(thirCol)
            fourCol = float(columns[3])
            lc_Pz.append(fourCol)
            fifCol = float(columns[4])
            lc_Pz_dev.append(fifCol)

wat_Pz = np.array(wat_Pz)
lc_Pz = np.array(lc_Pz)
wat_Pz = wat_Pz * 10000
lc_Pz = lc_Pz * 10000

total_Pz = wat_Pz + lc_Pz

plt.rcParams['axes.linewidth'] = 1.5
# plt.xticks(rotation=25)
plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.ylabel(r'$\mathrm{P_z(z)(10^{-4}\/e/\AA^2)}$', fontsize=27)
plt.tick_params(which='both', width=2)
plt.xlim((binCenter[0],binCenter[len(binCenter)-1]))
plt.ylim((-18, 18))
# plt.yticks(weight='bold')
# plt.xticks(weight='bold')
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(10, 5)
plt.plot(binCenter, wat_Pz, color='b', linewidth=2, label='Water')
plt.plot(binCenter, lc_Pz, color='r', linewidth=2, label='LC')
plt.plot(binCenter, total_Pz, color='g', linewidth=2, linestyle='dashed', label='Sum')
plt.legend(loc='best', frameon=True, fancybox=True)
plt.savefig('../Analysis/plot_Pz.png', dpi=600, facecolor='w', edgecolor='w',
            orientation='portrait', papertype='letter', format='png',
            transparent=True, bbox_inches='tight', pad_inches=0.1,
            frameon=None)
plt.clf()
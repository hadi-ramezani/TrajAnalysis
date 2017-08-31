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
                    action="store",nargs='?',
                    required=True, dest="dz", type=float,
                    help="Grid thickness")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")                                        
parser.add_argument("-o",
                    action="store", nargs='?', 
                    required=True, dest="output1",
                    help="output filename for ScalarZ")
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
                    action="store", nargs='?', 
                    required=True, dest="select",
                    help="atoms selection argument")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
dz = args.dz
output_filename = args.output1
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection = args.select


time1 = time.time()
output_file = open(output_filename,"w")

# Read the trj and psf file
if traj_filename == "" :
    list = [str(firstFile)+".dcd"]
    for i in xrange(firstFile+1,lastFile+1) :
        list.append(str(i)+".dcd")
    print "List of Trajectory files:" + str(list)
    u = Universe(psf_filename, list)
else :
    u = Universe(psf_filename, traj_filename)

# Obtain initial information form psf and trj files
natoms = len(u.selectAtoms(atom_selection))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)

# Define Slabs
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
  end_frame = num_frames
print "reading frame " + str(start_frame) + " to " + str(end_frame)

# Deine arrays

dens_f = np.empty((0),dtype='float')

# PBC
box = u.dimensions[0:3]
volume = dz * box[0] * box[1]

#Conversion coefficient (from a.u. to e/A^3)
conv_coeff = 1.67377

pvector = np.zeros((natoms),dtype='float')
nframes = 0
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
        coor = u.selectAtoms(atom_selection).coordinates()
        mass = u.selectAtoms(atom_selection).masses()
        intMass2 = mass * 0.5
        intMass2 = intMass2.astype(int)
        charge = u.selectAtoms(atom_selection).charges()
        massMinusCharge = intMass2 - charge


        #Calculate the position of each slab to have almost same number of atoms in each slab
        if curr_frame == start_frame :
            for i in xrange(0,len(coor)) :
                loc = coor[i]
                pvector[i] = loc[2]
            
            pvector = np.sort(pvector)
            highZ = int(pvector[natoms-1]) + 5
            lowZ = int(pvector[0])
            nslab = int((highZ - lowZ)/dz)
            bins = np.zeros((nslab+1),dtype='float')
            bins[nslab] = highZ
            bins[0] = lowZ
            slab_ind = nslab - 1
            dens = np.zeros((nslab),dtype='int')
            for i in xrange(1,nslab):
                bins[i] = bins[i-1] + dz

        dens.fill(0)
        for i in xrange(0,len(coor)) :
            loc = coor[i]
            p = loc[2]
            x = min(int((p-bins[0])/dz),slab_ind)
            dens[x] += massMinusCharge[i]

        dens_f = np.append(dens_f,dens,0)



dens_f = np.reshape(dens_f,( len(dens_f)/nslab,nslab))
dens_average =  np.average(dens_f,0)
dens_deviation =  np.std(dens_f,0)       
dens_average = conv_coeff * dens_average/volume
dens_deviation =  conv_coeff * dens_deviation/volume
p0 = (bins[0]+bins[1])/2
output_file.write('%5.3f %10.7f %10.7f\n' % ( bins[0] ,0, 0))
for x in range(nslab):
    p = (bins[x]+bins[x+1])/2
    d = (bins[x+1]-bins[x])
    output_file.write('%5.3f %10.7f %10.7f\n' % ( p , dens_average[x], dens_deviation[x]))
output_file.write('%5.3f %10.7f %10.7f\n' % ( bins[nslab] ,0, 0))


time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()

#Plot the electron density profiles
# Read the experimental results
outFile_exptl = open('../Analysis/electron_dens_exptl.dat', 'r')
binCenter_exptl=[]
density_exptl = []

for line in outFile_exptl:
    line = line.strip()
    if line != "":
        columns = line.split()
        if '#' in columns[0]:
            check='yes'
        else:
            firstCol = float(columns[0])
            binCenter_exptl.append(firstCol)
            secCol = float(columns[1])
            density_exptl.append(secCol)

# Read the simulated values
outFile = open(str(output_filename), 'r')
binCenter_sim=[]
density_sim = []
stDev_sim = []

for line in outFile:
    line = line.strip()
    if line != "":
        columns = line.split()
        if '#' in columns[0]:
            check='yes'
        else:
            firstCol = float(columns[0])
            binCenter_sim.append(firstCol)
            secCol = float(columns[1])
            density_sim.append(secCol)
            thirCol = float(columns[2])
            stDev_sim.append(thirCol)

# Convert all lists to numpy arrays
binCenter_exptl = np.array(binCenter_exptl); density_exptl = np.array(density_exptl)
binCenter_sim   = np.array(binCenter_sim); density_sim = np.array(density_sim); stDev_sim = np.array(stDev_sim)

# Select the members in the range z<150 and assign it to the same variables
binCenter_exptl = binCenter_exptl[np.where(binCenter_exptl < 150)]
density_exptl   = density_exptl[:len(binCenter_exptl)]

# Choose the last member of the experimental data (in the range z < 150) and call it bulk value
bulk_density_exptl = density_exptl[-1]
bin_bulk_exptl = binCenter_exptl[-1]

# Select the members in the range z<200 and assign it to the same variables (this is to remove the possibility of the max electron density near the other interface)
binCenter_sim = binCenter_sim[np.where(binCenter_sim < 200)]
density_sim   = density_sim[:len(binCenter_sim)]
stDev_sim     = stDev_sim[:len(binCenter_sim)]

# Find the maximum values of experimental and simulated electron density and the corresponding bin
max_density_exptl = np.max(density_exptl)
bin_max_exptl = binCenter_exptl[np.argmax(density_exptl)]

max_density_sim = np.max(density_sim)
bin_max_sim = binCenter_sim[np.argmax(density_sim)]

#Shift the binCenter_sim to match the max electron density of simulation and experiment
binCenter_sim = binCenter_sim - (bin_max_sim        - bin_max_exptl)

# Select the members in the range z<150 and assign it to the same variables
binCenter_sim = binCenter_sim[np.where(binCenter_sim < 150)]
density_sim   = density_sim[:len(binCenter_sim)]
stDev_sim     = stDev_sim[:len(binCenter_sim)]

# Choose the last member of the simulated data (in the range z < 150) and call it bulk value
bulk_density_sim = np.average(density_sim[-20:])
bin_bulk_sim = binCenter_sim[-1]

#Scale the simulation data in a way that it goes from zero at vacuum interface to 1 in the bulk
scale = 1/bulk_density_sim
density_sim   = density_sim * scale

#Plot the electron density data
plt.rcParams['axes.linewidth'] = 1.5
#plt.xticks(rotation=25)
plt.ylabel(r'$\mathrm{Electron\/density\/(a.u.)}$', fontsize=27)
plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.tick_params(which='both', width=2)
plt.ylim((0,2))
plt.xlim(-20,150)
#plt.yticks(weight='bold')
#plt.xticks(weight='bold')
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(8,5)
plt.errorbar(binCenter_sim, density_sim, xerr=stDev_sim, marker='o', markeredgecolor='b', label='Simulation')
plt.plot(binCenter_exptl, density_exptl, marker='o', markeredgecolor='g', label='Experiment')
plt.legend(loc='best', frameon=True, fancybox=True)
plt.savefig('../Analysis/plot_electron_density.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

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
parser.add_argument("-n",
                    action="store",nargs='?',
                    required=True, dest="nslab", type=int,
                    help="Number of slabs")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")                                        
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/mix_density.dat",
                    required=False, dest="output1",
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
parser.add_argument("-sn",
                    action="store", nargs='?', default="name NY1",
                    required=False, dest="atom_selection_N",
                    help="Atom selection argument for N")
parser.add_argument("-sc",
                    action="store", nargs='?', default="name CA12",
                    required=False, dest="atom_selection_C",
                    help="Atom selection argument for C")
parser.add_argument("-p5",
                    action="store", nargs='?', default=float('50.0'),
                    required=False, dest="initPerc5CB", type=float,
                    help="Mole percentage of 5CB")
parser.add_argument("-p8",
                    action="store", nargs='?', default=float('50.0'),
                    required=False, dest="initPerc8CB", type=float,
                    help="Mole percentage of 8CB")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
nslab = args.nslab
output_filename1 = args.output1
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
initPerc5CB = args.initPerc5CB
initPerc8CB = args.initPerc8CB

atom_selection_N = args.atom_selection_N
atom_selection_C = args.atom_selection_C
atom_selection_N5 = "resname 5CB and name NY1"
atom_selection_N8 = "resname 8CB and name NY1"

time1 = time.time()

def position(v1,v2):
    return (v1+v2)/2.0

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
natoms = len(u.selectAtoms(atom_selection_N))
natoms5 = len(u.selectAtoms(atom_selection_N5))
natoms8 = len(u.selectAtoms(atom_selection_N8))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)
print "Number of 5CB molecules = " + str(natoms5)
print "Number of 8CB molecules = " + str(natoms8)

# Define Slabs
delta = int(natoms/nslab)
print "Number of atoms in each slab is around " + str(delta)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
  end_frame = num_frames
print "reading frame " + str(start_frame) + " to " + str(end_frame)

# Deine arrays
bins = np.zeros((nslab+1),dtype='float')

pvector = np.zeros((natoms),dtype='float')

dens = np.zeros((nslab), dtype='float')
dens5 = np.zeros((nslab), dtype='float')
dens8 = np.zeros((nslab), dtype='float')
dens_f = np.empty((0), dtype='float')
dens5_f = np.empty((0), dtype='float')
dens8_f = np.empty((0), dtype='float')

nframes = 0
# Loop over snapshots and calculate density
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
        coor_N = u.selectAtoms(atom_selection_N).coordinates()
        coor_C = u.selectAtoms(atom_selection_C).coordinates()

        coor_N5 = u.selectAtoms(atom_selection_N5).coordinates()
        coor_N8 = u.selectAtoms(atom_selection_N8).coordinates()

        if len(coor_N) != len(coor_C) :
          print >>sys.stderr, "Error: number of atoms in each group does not match"
          sys.exit(1)
        
        #Calculate the position of each slab to have almost same number of atoms in each slab
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            pvector[i] = position(loc_C[2],loc_N[2])
            
        pvector = np.sort(pvector)
        bins[nslab] = pvector[natoms-1]
        for i in range(nslab):
            bins[i] = pvector[i*delta]
        
        dens.fill(0)
        dens5.fill(0)
        dens8.fill(0)
        j_5CB = 0
        j_8CB = 0
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            p = position(loc_C[2],loc_N[2])
            x = 1
            # Due to dynamics binning it is not possible to get
            # the bin index directly from z coordinate, hence using while loop
            while x<=nslab :
              if p<=bins[x] :
                 dens[x-1] += 1
                 break
              else :
                 x += 1
            if (j_5CB < natoms5 and np.allclose(loc_N, coor_N5[j_5CB])) :
                j_5CB += 1
                dens5[x-1] += 1
            if (j_8CB < natoms8 and np.allclose(loc_N, coor_N8[j_8CB])) :
                j_8CB += 1
                dens8[x-1] += 1
        dens_f = np.append(dens_f, dens, 0)
        dens5_f = np.append(dens5_f, dens5, 0)
        dens8_f = np.append(dens8_f, dens8, 0)
        
dens_f = np.reshape(dens_f,(len(dens_f)/nslab,nslab))
dens_average =  np.average(dens_f,0)
dens_deviation =  np.std(dens_f,0)

dens5_f = np.reshape(dens5_f,(len(dens5_f)/nslab,nslab))
dens5_average =  np.average(dens5_f,0)
dens5_deviation =  np.std(dens5_f,0)

dens8_f = np.reshape(dens8_f,(len(dens8_f)/nslab,nslab))
dens8_average =  np.average(dens8_f,0)
dens8_deviation =  np.std(dens8_f,0)

# Write the data
p = (bins[:-1] + bins[1:])/2
np.savetxt(output_filename1,np.c_[p, dens5_average, dens5_deviation,
                                  dens8_average, dens8_deviation], fmt='%5.3f %5.3f %5.3f %5.3f %5.3f')

time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)

#Read-in the data
outFile = open(str(output_filename1), 'r')
binCenter = []
dens5_avg = []; dens5_std = []
dens8_avg = []; dens8_std = []

for line in outFile:
    line = line.strip()
    if line != "":
        columns = line.split()
        if not line.startswith('#'):
            binCenter.append(float(columns[0]))
            dens5_avg.append(float(columns[1]))
            dens5_std.append(float(columns[2]))
            dens8_avg.append(float(columns[3]))
            dens8_std.append(float(columns[4]))

# Convert the lists to numpy arrays
dens5_avg = np.array(dens5_avg); dens5_std = np.array(dens5_std)
dens8_avg = np.array(dens8_avg); dens8_std = np.array(dens8_std)

# Convert the numbers to percentages
perc_5cb = 100 * dens5_avg/(dens5_avg + dens8_avg)
std_5cb = 100 * dens5_std/(dens5_avg + dens8_avg)
perc_8cb = 100 * dens8_avg/(dens5_avg + dens8_avg)
std_8cb = 100 * dens8_std/(dens5_avg + dens8_avg)

#Plot the density data
plt.rcParams['axes.linewidth'] = 1.5
plt.xticks(rotation=25)	
plt.xlabel(r'$\mathrm{Composition}\/(\%)$', fontsize=27, weight='bold')
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)	
plt.tick_params(which='both', width=2)	
plt.xlim((0,100))
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(5,10)
plt.plot(perc_5cb, binCenter, marker='o', markersize = 4, color = 'b', markeredgecolor = 'b', label = '5CB')
plt.plot(perc_8cb, binCenter, marker='o', markersize = 4, color = 'r', markeredgecolor='r', label = '8CB')
plt.fill_betweenx(binCenter, perc_5cb + std_5cb, perc_5cb - std_5cb,
                 alpha = 0.2, facecolor = 'b', edgecolor = 'none')
plt.fill_betweenx(binCenter, perc_8cb + std_8cb, perc_8cb - std_8cb,
                 alpha = 0.2, facecolor = 'r', edgecolor = 'none')
plt.axvline(x = initPerc5CB, linestyle='--', linewidth=2, color='b')
plt.axvline(x = initPerc8CB, linestyle='--', linewidth=2, color='r')
plt.legend(loc='center left', frameon=True, fancybox=True, fontsize=20)

plt.savefig('../Analysis/mix_dens.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

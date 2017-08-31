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
                    action="store", nargs='?', default="../Analysis/ScalerZ.dat", 
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
nslab = args.nslab
output_filename = args.output
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName

atom_selection_N = "name NY1"
atom_selection_C = "name CA12"
max_distance = 15 #The distance between two ends of the molecular vector must not be large than this value

time1 = time.time()
output_file = open(output_filename,"w")
#output_file2 = open("Pz.dat","w")

def distance(v1, v2) :
    return sqrt((v1[0] - v2[0])*(v1[0] - v2[0]) + (v1[1] - v2[1])*(v1[1] - v2[1]) + (v1[2] - v2[2])*(v1[2] - v2[2]) )
def distance_d(v1, v2) :
    return abs(v1-v2)
def vector(v1, v2):
    return [v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]
def position(v1,v2):
    return (v1+v2)/2.0
def Sz(v) : 
    return (3.*v[2]*v[2]-1.)/2.0
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
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)

# Define Slabs
resolution = 10
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
Scaler_Z = np.empty((nslab),dtype='float')
Scaler_Z_f = np.empty((0),dtype='float')

dens = np.zeros((nslab),dtype='int')
dens_p = np.zeros((nslab,resolution+1), dtype='int')

Pro = np.zeros((nslab,resolution+1), dtype='float')
Pro_t = np.zeros((nslab,resolution+1), dtype='float')

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
        coor_N = u.selectAtoms(atom_selection_N).coordinates()
        coor_C = u.selectAtoms(atom_selection_C).coordinates()
#        print np.amin(coor_N)
        if len(coor_N) != len(coor_C) :
          print >>sys.stderr, "Error: number of atoms in each group does not match"
          sys.exit(1)
        
        # Applay PBC 
        box = u.dimensions[0:3]
#        box2 = box/2
#        for i in xrange(0,len(coor_N)) :
#            loc_N = coor_N[i]
#            loc_C = coor_C[i]
#            dis = distance(loc_N, loc_C)
#            if dis > max_distance :
#               for j in xrange(0,len(loc_N)):
#                  dis_d = distance_d(loc_N[j],loc_C[j])
#                  if dis_d >= box2[j] :
#                     if loc_N[j] > box2[j] :
#                        loc_N[j] = loc_N[j] - box[j]
#                     else :
#                        loc_N[j] = loc_N[j] + box[j]
#                  dis_d = distance_d(loc_N[j],loc_C[j])
#            dis = distance(loc_N, loc_C)
#            if dis > max_distance :
#               print >>sys.stderr, "Error"
        
        #Calculate the opsition of each slab to have almost same number of atoms in each slab
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            pvector[i] = position(loc_C[2],loc_N[2])
            
        pvector = np.sort(pvector)
        bins[nslab] = pvector[natoms-1]
        for i in range(nslab):
            bins[i] = pvector[i*delta]
        
        dens.fill(0)
        Scaler_Z.fill(0)
        dens_p.fill(0)
        Pro.fill(0)
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            p = position(loc_C[2],loc_N[2])
            x = 1
            while x<=nslab :
              if p<=bins[x] :
                 vv = vector(loc_C,loc_N)
                 vv = vv/np.linalg.norm(vv)
                 Scaler = Sz(vv)
                 Scaler_Z[x-1] += Scaler
                 dens[x-1] += 1
                 y = int(np.absolute(Scaler)*resolution)
                 dens_p[x-1,y] += 1
                 Pro[x-1,y] += 1

                 break
              else :
                 x += 1
        Scaler_Z = Scaler_Z / dens
        Scaler_Z_f = np.append(Scaler_Z_f,Scaler_Z,0)
        
        for xx in range(nslab): 
            for yy in range(resolution+1):
                Pro_t[xx,yy] += Pro[xx,yy] / natoms

Scaler_Z_f = np.reshape(Scaler_Z_f,( len(Scaler_Z_f)/nslab,nslab))
S_average =  np.average(Scaler_Z_f,0)
S_deviation =  np.std(Scaler_Z_f,0)
p0 = (bins[0]+bins[1])/2
for x in range(nslab):
    p = (bins[x]+bins[x+1])/2
    d = (bins[x+1]-bins[x])
    output_file.write('%5.3f %5.3f %5.3f %10.7f\n' % ( p ,S_average[x], S_deviation[x], dens[x]/box[0]/box[1]/d))

time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()

#Plot the output results
outFile = open(str(output_filename), 'r')
binCenter=[]
orderParam = []
stdDev = []
density = []
for line in outFile:
	line = line.strip()
	if line != "":
		columns = line.split()
		if '#' in columns[0]:
			check='yes'
		else:
			firstCol = float(columns[0])
			binCenter.append(firstCol)
			secCol = float(columns[1])
			orderParam.append(secCol) 
			thirCol = float(columns[2])
			stdDev.append(thirCol)
			fourCol = float(columns[3])
			density.append(fourCol)

#Plot the P2 data
x=np.zeros(len(binCenter))
plt.plot(x, binCenter, '--', linewidth=1, color='red')				
plt.xlabel(r'$\mathrm{P_2}$', fontsize=16)
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=16)	
plt.tick_params(which='both', width=2)		
plt.xlim((-0.5,1))
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
#plt.yticks(weight='bold')
#plt.xticks(weight='bold')
plt.gcf().set_size_inches(5,10)
plt.errorbar(orderParam, binCenter, xerr=stdDev, marker='o', markeredgecolor='b')
plt.savefig('../Analysis/plot_P2.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
plt.clf()		
#Plot the density data	
plt.xticks(rotation=25)	
plt.xlabel(r'$\mathrm{density\/(molecule/} \AA^3)$', fontsize=16, weight='bold')
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=16)	
plt.tick_params(which='both', width=2)	
gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%.e'))	
#plt.xlim((-0.5,1))
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
#plt.yticks(weight='bold')
#plt.xticks(weight='bold')
plt.gcf().set_size_inches(5,10)
plt.plot(density, binCenter, marker='o', markeredgecolor='b')
plt.savefig('../Analysis/plot_dens.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

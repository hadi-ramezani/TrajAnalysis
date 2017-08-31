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
                    action="store", nargs='?', default="../Analysis/tensorZ.dat", 
                    required=False, dest="output",
                    help="output filename")
parser.add_argument("-b",
                    action="store", nargs='?', default=int('0'), 
                    required=False, dest="start_frame", type=int,
                    help="First frame to start reading trajectory")
parser.add_argument("-e",
                    action="store", nargs='?', default=float('inf'), 
                    required=False, dest="end_frame", type=int,
                    help="Last frame to stop readining trajectory")
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


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
nslab = args.nslab
output_filename = args.output
start_frame = args.start_frame
end_frame = args.end_frame
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection_N = args.atom_selection_N
atom_selection_C = args.atom_selection_C

max_distance = 15 #The distance between two ends of the molecular vector must not be large than this value

time1 = time.time()
output_file = open(output_filename + "_local","w")
output_file2 = open(output_filename + "_global","w")
output_file3 = open(output_filename + "_vector","w")

def distance(v1, v2) :
    return sqrt((v1[0] - v2[0])*(v1[0] - v2[0]) + (v1[1] - v2[1])*(v1[1] - v2[1]) + (v1[2] - v2[2])*(v1[2] - v2[2]) )
def distance_d(v1, v2) :
    return abs(v1-v2)
def vector(v1, v2):
    return [[v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]]
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
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)

# Define Slabs
delta = int(natoms/nslab)
print "Number of atoms in each slab is around " + str(delta)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
   end_frame = num_frames
print "read frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays
bins = np.zeros((nslab+1),dtype='float')

pvector = np.zeros((natoms),dtype='float')
Q_slab = np.empty((nslab,9),dtype='float')
Qt = np.zeros((3,3),dtype='float')
Scaler = np.empty((nslab),dtype='float')
dataout = np.zeros((5),dtype='float')
dataout_t = np.empty((0),dtype='float')
Scaler_f = np.empty((0),dtype='float')
vector_x = np.empty((nslab),dtype='float')
vector_x_f = np.empty((0),dtype='float')
vector_y = np.empty((nslab),dtype='float')
vector_y_f = np.empty((0),dtype='float')
vector_z = np.empty((nslab),dtype='float')
vector_z_f = np.empty((0),dtype='float')
dens = np.zeros((nslab),dtype='int')
I = np.matrix([[1.,0,0],[0,1.,0],[0,0,1.]])

nframes = 0
# Loops in the trj
for curr_frame in xrange(0, num_frames) :
    if curr_frame != 0 :
        trj = u.trajectory.next()
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame >= start_frame and curr_frame <= end_frame :
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
        coor_N = u.selectAtoms(atom_selection_N).coordinates()[:natoms]
        coor_C = u.selectAtoms(atom_selection_C).coordinates()[:natoms]
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
        # Calculate Q tensor for atoms
        # Calculate Q tensor for each slab by summing Q of atomes located in the slab
        # Calculate global Q 
        # Calculate number of atoms in each slab
        dens.fill(0)
        Q_slab.fill(0)
        Qt.fill(0)
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            p = position(loc_C[2],loc_N[2])
            x = 1
            while x<=nslab :
              if p<=bins[x] :
                 vv = vector(loc_C,loc_N)
                 vv = vv/np.linalg.norm(vv)
                 Q = (3*np.dot(vv.T,vv) - I)/2.0 # factor 3/2 means that biggest eigvalu = nematic order
                 Qt += Q
                 Q = np.reshape(Q,(9))
                 Q_slab[x-1,:] = Q + Q_slab[x-1,:]
                 dens[x-1] += 1
                 break
              else :
                 x += 1
        for i in range(nslab):
            Q = Q_slab[i,:] / dens[i] # Normalize local Q
            Q = np.reshape(Q,(3,3))
            e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
#            print e_value,  np.argmax(e_value), np.amax(e_value)
            Scaler[i] = np.amax(e_value)   # Larger eigenvalue is S
            vector_x[i] = np.absolute(e_vector[0,np.argmax(e_value)])
            vector_y[i] = np.absolute(e_vector[1,np.argmax(e_value)])
            vector_z[i] = np.absolute(e_vector[2,np.argmax(e_value)])
            
        Qt = Qt/natoms  # Normalize global Q
        e_value, e_vector = np.linalg.eigh(Qt)  # Calculate global order parameter
        # Store time, Max_eigenvalue and eigenvector
        dataout[0] = curr_time
        dataout[1] = np.amax(e_value)
        dataout[2:5] = e_vector[:,np.argmax(e_value)]
        dataout_t = np.append(dataout_t, dataout, 0)
        # Store local S and eigenvector
        Scaler_f = np.append(Scaler_f, Scaler,0)
        vector_x_f = np.append(vector_x_f, vector_x, 0)
        vector_y_f = np.append(vector_y_f, vector_y, 0)
        vector_z_f = np.append(vector_z_f, vector_z, 0)

#Write output file
dataout_t = np.reshape(dataout_t,(len(dataout_t)/5,5))
print >> output_file2, "@ time(ps)  Max_eign  eigen_vector "
np.savetxt(output_file2, dataout_t,fmt='%7.4f')
ss = np.average(dataout_t,0)
st = np.std(dataout_t,0)
print " Mean value of global order parameter = " + str(ss[1]) + " STD " + str(st[1])

# Calculate mean and STD for local S
Scaler_f = np.reshape(Scaler_f,( len(Scaler_f)/nslab,nslab))
S_average =  np.average(Scaler_f,0)
S_deviation =  np.std(Scaler_f,0)
# Calculate mean value of eigen vectro vector
vector_x_f = np.reshape(vector_x_f,(len(vector_x_f)/nslab,nslab))
vector_x_av = np.average(vector_x_f,0)
vector_x_dev = np.std(vector_x_f,0)

vector_y_f = np.reshape(vector_y_f,(len(vector_y_f)/nslab,nslab))
vector_y_av = np.average(vector_y_f,0)
vector_y_dev = np.std(vector_y_f,0)

vector_z_f = np.reshape(vector_z_f,(len(vector_z_f)/nslab,nslab))
vector_z_av = np.average(vector_z_f,0)
vector_z_dev = np.std(vector_z_f,0)
#write output file
p0 = (bins[0]+bins[1])/2
for x in range(nslab):
    p = (bins[x]+bins[x+1])/2
    output_file.write('%5.3f %5.3f %5.3f \n' % ( p ,S_average[x] ,S_deviation[x]))
    output_file3.write('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n' % ( p , vector_x_av[x], vector_x_dev[x], vector_y_av[x], vector_y_dev[x], vector_z_av[x], vector_z_dev[x]))
time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()
output_file3.close()

#Plot the output results
outFile = open(str(output_filename + "_local"), 'r')
binCenter=[]
orderParam = []
stdDev = []
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

#Plot the S data
plt.rcParams['axes.linewidth'] = 1.5
plt.xticks(rotation=25)
plt.xlabel(r'$\mathrm{S}$', fontsize=27)
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)	
plt.tick_params(which='both', width=2)		
plt.xlim((0,1))
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
#plt.yticks(weight='bold')
#plt.xticks(weight='bold')
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(5,10)
plt.errorbar(orderParam, binCenter, xerr=stdDev, color='b', marker='o', markeredgecolor='b')
plt.savefig('../Analysis/plot_S.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
plt.clf()

#Plot the director componenets
outFile = open(str(output_filename + "_vector"), 'r')
binCenter=[]
vx = []
vx_dev =[]
vy = []
vy_dev = []
vz = []
vz_dev = []
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
			vx.append(secCol) 
			thirCol = float(columns[2])
			vx_dev.append(thirCol)
			fourCol = float(columns[3])
			vy.append(fourCol)
			fifCol = float(columns[4])
			vy_dev.append(fifCol)
			sixCol = float(columns[5])
			vz.append(sixCol)
			sevenCol = float(columns[6])
			vz_dev.append(sevenCol)
						
#Plot the nk data
plt.rcParams['axes.linewidth'] = 1.5
plt.xticks(rotation=25)
plt.xlabel(r'$\mathrm{n_k}$', fontsize=27)
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)	
plt.tick_params(which='both', width=2)		
plt.xlim((0,1))
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
#plt.yticks(weight='bold')
#plt.xticks(weight='bold')
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(5,10)
plt.errorbar(vx, binCenter, xerr=vx_dev, color='r', marker='o', markeredgecolor='r', label=r'$\mathrm{n_x}$')
plt.errorbar(vy, binCenter, xerr=vy_dev, color='g', marker='o', markeredgecolor='g', label=r'$\mathrm{n_y}$')
plt.errorbar(vz, binCenter, xerr=vz_dev, color='m', marker='o', markeredgecolor='m', label=r'$\mathrm{n_z}$')
plt.legend(loc='best', frameon=True, fancybox=True)
plt.savefig('../Analysis/plot_nk.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

		

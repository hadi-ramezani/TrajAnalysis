#! /usr/bin/env python

from MDAnalysis import *
from math import *
import numpy as np
from scipy import interpolate
import time
import sys
import argparse
import os

parser = argparse.ArgumentParser(description="Read trajectory and coordinate file")

parser.add_argument("-t", 
                    action="store", nargs='?', default="",
                    required=False, dest="traj", 
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-nx",
                    action="store",nargs='?',
                    required=True, dest="nslabx", type=int,
                    help="Number of slabs")
parser.add_argument("-ny",
                    action="store",nargs='?',
                    required=True, dest="nslaby", type=int,
                    help="Number of slabs")
parser.add_argument("-vx",
                    action="store",nargs='?',
                    required=True, dest="vslabx", type=int,
                    help="Number of slabs")
parser.add_argument("-vy",
                    action="store",nargs='?',
                    required=True, dest="vslaby", type=int,
                    help="Number of slabs")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/Output", 
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
nslabx = args.nslabx
nslaby = args.nslaby
vslabx = args.vslabx
vslaby = args.vslaby
output_filename = args.output
start_frame = args.start_frame
end_frame = args.end_frame
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection_N = args.atom_selection_N
atom_selection_C = args.atom_selection_C


time1 = time.time()
output_file = open(output_filename + "_scaler","w")
output_file2 = open(output_filename + "_vector" ,"w")

def vector(v1, v2):
    return [[v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]]
def position(v1,v2):
    return (v1+v2)/2.0
def anint(a):
    if a > 0 :
       an = int(a+0.5)
    else :
       an = int(a-0.5)
    return int(an)
# Read the trj and psf file
if traj_filename == "" :
  list = [str(firstFile)+".dcd"]
  for i in xrange(firstFile+1,lastFile+1) :
    list.append(str(i)+".dcd")
  print "List of Trajectory files:" + str(list)
  u = Universe(psf_filename, list)
else :
  u = Universe(psf_filename, traj_filename)

# Obtain initial information form gro and trj files
natoms = len(u.selectAtoms(atom_selection_N))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)
delta = (natoms/nslabx/nslaby)
deltav = (natoms/vslabx/vslaby)
print "Number of atoms in each slab is around " + str(delta)
print "Number of atoms in each slab is around " + str(deltav)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
   end_frame = num_frames
print "read frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays
binsx = np.zeros((nslabx+1),dtype='float')
binsy = np.zeros((nslaby+1),dtype='float')
vbinsx = np.zeros((vslabx+1),dtype='float')
vbinsy = np.zeros((vslaby+1),dtype='float')

Q_slab = np.empty((nslabx,nslaby,9),dtype='float')
dens = np.zeros((nslabx,nslaby),dtype='int')
dens_f = np.zeros((nslabx,nslaby),dtype='float')
Scalar = np.empty((nslabx,nslaby),dtype='float')

Q_slabv = np.empty((vslabx,vslaby,9),dtype='float')
densv = np.zeros((vslabx,vslaby),dtype='int')
densv_f = np.zeros((vslabx,vslaby),dtype='int')
Scalarv = np.empty((vslabx,vslaby),dtype='float')

vector_x = np.zeros((vslabx,vslaby),dtype='float')
vector_y = np.zeros((vslabx,vslaby),dtype='float')
vector_z = np.zeros((vslabx,vslaby),dtype='float')

Qt = np.zeros((3,3),dtype='float')
I = np.matrix([[1.,0,0],[0,1.,0],[0,0,1.]])

Scalar.fill(0)
Scalarv.fill(0)
Q_slab.fill(0)
Q_slabv.fill(0)
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
        coor_N = u.selectAtoms(atom_selection_N).coordinates()
        coor_C = u.selectAtoms(atom_selection_C).coordinates()
        if len(coor_N) != len(coor_C) :
           print >>sys.stderr, "Error: number of atoms in each group does not match"
           sys.exit(1)
        
        #Calculate the opsition of each slab to have almost same number of atoms in each slab
        #The method is correct if the falctuation of volume is very small
        box = u.dimensions[0:3]
        box2 = box/2
        dx = box[0]/nslabx
        dy = box[1]/nslaby
        volume = (dx * dy * box[2]) / 1000  # nm^3
        for i in range(nslabx+1):
            binsx[i] = i*dx - box2[0]
        for i in range(nslaby+1):
            binsy[i] = i*dy - box2[1]
 
        dxv = box[0]/vslabx
        dyv = box[1]/vslaby
        volumev = (dxv * dyv * box[2]) / 1000  # nm^3
        for i in range(vslabx+1):
            vbinsx[i] = i*dxv - box2[0]
        for i in range(vslaby+1):
            vbinsy[i] = i*dyv - box2[1]
       
#        # Calculate Q tensor for atoms
#        # Calculate Q tensor for each slab by summing Q of atomes located in the slab
#        # Calculate global Q 
#        # Calculate number of atoms in each slab
        Qt.fill(0)
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]

            px = position(loc_C[0],loc_N[0])
            py = position(loc_C[1],loc_N[1])

            # apply periodic boundary conditions
            px = px - box[0]*anint(px/box[0]) + box2[0]
            py = py - box[1]*anint(py/box[1]) + box2[1]
            
            vv = vector(loc_C,loc_N)
            vv = vv/np.linalg.norm(vv)
            Q = (3*np.dot(vv.T,vv) - I)/2.0 # factor 3/2 means that biggest eigvalu = nematic order
            Qt += Q
            Q = np.reshape(Q,(9))
           
            nx = int(px/dx) 
            ny = int(py/dy)
            if nx == nslabx : nx -= 1
            if nx < 0 :       nx += 1
            if ny == nslaby : ny -= 1
            if ny < 0 :       ny += 1
            dens[nx,ny] += 1
            Q_slab[nx,ny,:] = Q + Q_slab[nx,ny,:]

            nx = int(px/dxv) 
            ny = int(py/dyv)
            if nx == vslabx : nx -= 1
            if nx < 0 :       nx += 1
            if ny == vslaby : ny -= 1
            if ny < 0 :       ny += 1
            densv[nx,ny] += 1
            Q_slabv[nx,ny,:] = Q + Q_slabv[nx,ny,:]

        Qt = Qt/natoms  # Normalize global Q
        e_value, e_vector = np.linalg.eigh(Qt)  # Calculate global order parameter
        print "time " + str(curr_time) + " Scalar " + str(np.amax(e_value))

print "Total number of frames " + str(nframes)

for x in range(nslabx):
  for y in range(nslaby):
      dens_f[x,y] = float(dens[x,y]) / nframes
      if dens_f[x,y]>delta*0.2 :
         Q = Q_slab[x,y,:] / dens[x,y]  
         Q = np.reshape(Q,(3,3))
         e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
         Scalar[x,y] = np.amax(e_value)   # Larger eigenvalue is S
      else : 
         Scalar[x,y] = 0.6

for x in range(vslabx):
  for y in range(vslaby):
      densv_f[x,y] = float(densv[x,y]) / nframes
      if densv_f[x,y]>deltav*0.2 :
         Q = Q_slabv[x,y,:] / densv[x,y]  
         Q = np.reshape(Q,(3,3))
         e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
         Scalarv[x,y] = np.amax(e_value)   # Larger eigenvalue is S
         vector_x[x,y] = e_vector[0,np.argmax(e_value)]
         vector_y[x,y] = e_vector[1,np.argmax(e_value)]
         vector_z[x,y] = e_vector[2,np.argmax(e_value)]
      else : 
         Scalarv[x,y] = 0.6

##write Scalar
for x in range(nslabx):
  for y in range(nslaby):
    output_file.write('%5.3f %5.3f %5.3f %5.3f %5.3f \n' % ((binsx[x]+binsx[x+1])/2,(binsy[y]+binsy[y+1])/2,Scalar[x,y],dens_f[x,y],dens_f[x,y]/volume))
  output_file.write('\n')
##write vectors
zz = 0.0
for x in range(vslabx):
  for y in range(vslaby):
      output_file2.write('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n' % (((vbinsx[x+1]+vbinsx[x])/2.0),((vbinsy[y+1]+vbinsy[y])/2.0),zz,Scalarv[x,y],densv_f[x,y],vector_x[x,y],vector_y[x,y],vector_z[x,y]))
  output_file2.write('\n')
#
#
time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()
output_file2.close()
#

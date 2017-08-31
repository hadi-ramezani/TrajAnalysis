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
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/dwf.dat", 
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
parser.add_argument("-s",
                    action="store", nargs='?', default="circle",
                    required=False, dest="shape",
                    help="shape of the bubble: circle, cube, squircle")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
nslabx = args.nslabx
nslaby = args.nslaby
output_filename = args.output
start_frame = args.start_frame
end_frame = args.end_frame
firstFile = args.firstFileName
lastFile = args.lastFileName
shape = args.shape


atom_selection = "all"

time1 = time.time()
output_file = open(output_filename ,"w")

def anint(a):
    if a > 0 :
       an = int(a+0.5)
    else :
       an = int(a-0.5)
    return int(an)
def vector(v1, v2):
    return [[v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]]
def position(v1,v2):
    return (v1+v2)/2.0
def write_data():
    #write the data
    output_file.write('# x     y     3*dwf_x     3*dwf_y     3*dwf_z     dwf_sum     3*dwf_pr    1.5*dwf_or\n')
    for x in range(nslabx):
        for y in range(nslaby):
            output_file.write('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n' % (((binsx[x+1]+binsx[x])/2.0), ((binsy[y+1]+binsy[y])/2.0), 3*dwf_grid[x,y,0], 3*dwf_grid[x,y,1], 3*dwf_grid[x,y,2],
                                                                                   sum(dwf_grid[x,y,:]), 3*pr_dwf_grid[x,y], 1.5*or_dwf_grid[x,y]))
        output_file.write('\n')
def get_normal_direction(x, y):
    ''' This function calculate the director using gradient of the
     shape'''

    n = np.zeros(3)
    if shape == "circle":
        n[0] = x
        n[1] = y
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    elif shape == "cube":
        n[0] = np.power(x, 39)
        n[1] = np.power(y, 39)
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    elif shape == "squircle":
        n[0] = np.power(x, 3)
        n[1] = np.power(y, 3)
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    else:
        print sys.stderr,"Unknown shape, exiting ..."
        sys.exit(1)
    return n

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
natoms = len(u.selectAtoms(atom_selection))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)
delta = (natoms/nslabx/nslaby)
print "Number of atoms in each slab is around " + str(delta)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
   end_frame = num_frames
print "read frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays
binsx = np.zeros((nslabx+1),dtype='float')
binsy = np.zeros((nslaby+1),dtype='float')

dwf_grid = np.zeros((nslabx, nslaby,3),dtype='float')
pr_dwf_grid = np.zeros((nslabx, nslaby),dtype='float')
or_dwf_grid = np.zeros((nslabx, nslaby),dtype='float')
dens = np.zeros((nslabx, nslaby),dtype='int')
dens_f = np.zeros((nslabx, nslaby),dtype='float')


DWF = np.zeros(3, dtype='float')

box = u.dimensions[0:3]
box2 = box/2
# To make understanding of array indexes easier
x = 0; y = 1; z = 2
num_sele_frames = end_frame - start_frame

# A 3-dimensional array to store the wrapped coordinatescd
wrapped_coords = np.zeros((natoms,num_sele_frames,3),dtype='float')
dv_index = np.zeros((natoms,num_sele_frames,3),dtype='float')
unwrapped_coords = np.zeros((natoms,num_sele_frames,3),dtype='float')


# Go through trajectories and read wrapped trajectories
nframes = 0
for curr_frame in xrange(0, num_frames) :
    if curr_frame != 0 :
        trj = u.trajectory.next()
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame >= start_frame and curr_frame < end_frame :
        wrapped_coords[:,nframes,:] = u.selectAtoms(atom_selection).coordinates()
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
    elif curr_frame > end_frame:
        u.trajectory.close()
        break

# Unwrap Trajectories
print "Unwrapping trajectories"
dv0_wrapped     = np.zeros((natoms,1,3),dtype='float')
dv_rest_wrapped = wrapped_coords[:,1:,:] - wrapped_coords[:,:-1,:]
dv_wrapped      = np.concatenate((dv0_wrapped,dv_rest_wrapped),axis=1)

print "Deleting some of the unnecessary arrays to free memory"
del dv0_wrapped, dv_rest_wrapped

print "Computing indexes to unwrap the coordinates"
dv_index[:, :, x] = np.round(dv_wrapped[:, :, x] / box[x])
dv_index[:, :, y] = np.round(dv_wrapped[:, :, y] / box[y])
dv_index[:, :, z] = np.round(dv_wrapped[:, :, z] / box[z])


# One loop that uses the indexes to unwrap the coordinates
for frame in range(num_sele_frames):
    print "Unwrapping frame " + str(frame+1)
    unwrapped_coords[:,frame,x] = wrapped_coords[:,frame,x] - box[x] * np.sum(dv_index[:,:frame+1,x], axis=1)
    unwrapped_coords[:,frame,y] = wrapped_coords[:,frame,y] - box[y] * np.sum(dv_index[:,:frame+1,y], axis=1)
    unwrapped_coords[:,frame,z] = wrapped_coords[:,frame,z] - box[z] * np.sum(dv_index[:,:frame+1,z], axis=1)


print "Deleting some of the unnecessary arrays to free memory"
del dv_wrapped
dv0_unwrapped     = np.zeros((natoms,1,3),dtype='float')
dv_rest_unwrapped = unwrapped_coords[:,1:,:] - unwrapped_coords[:,:-1,:]
dv_unwrapped      = np.concatenate((dv0_unwrapped,dv_rest_unwrapped),axis=1)
dist = np.power(dv_unwrapped,2)
print "Deleting some other unnecessary arrays to free memory"
del dv0_unwrapped, dv_rest_unwrapped, unwrapped_coords

#Calculate the position of each slab to have almost same number of atoms in each slab
#The method is correct if the fluctuation of volume is very small
dx = box[0]/nslabx
dy = box[1]/nslaby
volume = (dx * dy * box[2]) / 1000  # nm^3
for i in range(nslabx+1):
    binsx[i] = i*dx - box2[0]
for i in range(nslaby+1):
    binsy[i] = i*dy - box2[1]

for frame in range(1, num_sele_frames, 1):

    DWF.fill(0)
    for i in xrange(0,len(wrapped_coords)) :
        loc = wrapped_coords[i,frame,:]

        px = loc[x]
        py = loc[y]

        n = get_normal_direction(px, py)

        # apply periodic boundary conditions
        px = px - box[0]*anint(px/box[0]) + box2[0]
        py = py - box[1]*anint(py/box[1]) + box2[1]

        dwf = dist[i,frame,:]

        pr_v = dv_unwrapped[i, frame, :] * n
        pr_dist = np.power(pr_v, 2)
        pr_dwf = np.sum(pr_dist)

        or_dist = dwf - pr_dist
        or_dwf = np.sum(or_dist)

        DWF += dwf
        nx = int(px/dx)
        ny = int(py/dy)

        if nx == nslabx : nx -= 1
        if nx < 0 :       nx += 1
        if ny == nslaby : ny -= 1
        if ny < 0 :       ny += 1

        dens[nx,ny] += 1
        dwf_grid[nx,ny,:] = dwf + dwf_grid[nx,ny,:]
        pr_dwf_grid[nx,ny] = pr_dwf + pr_dwf_grid[nx,ny]
        or_dwf_grid[nx,ny] = or_dwf + or_dwf_grid[nx,ny]

    DWF = DWF/natoms  # Normalize global DWF
    print " DWF " + str(np.linalg.norm(DWF))


for x in range(nslabx):
  for y in range(nslaby):
      dens_f[x,y] = float(dens[x,y]) / num_sele_frames
      if dens_f[x,y] > delta*0.2 :
          dwf_grid[x,y,:] = dwf_grid[x,y,:] / dens[x,y]
          pr_dwf_grid[x,y] = pr_dwf_grid[x,y] / dens[x,y]
          or_dwf_grid[x,y] = or_dwf_grid[x,y] / dens[x,y]
      else:
          dwf_grid[x,y,:] = 0
          pr_dwf_grid[x,y] = 0
          or_dwf_grid[x,y] = 0

write_data()
############################################################################################

time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()
#

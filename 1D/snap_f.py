#! /usr/bin/env python

from MDAnalysis import *
from math import *
import argparse
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Make triangle")

parser.add_argument("-t",
                    action="store", nargs='?',
                    required=True, dest="traj",
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-r",
                    action="store", nargs='?', default="head",
                    required=False, dest="head",
                    help="the upper bond of your scale bar")
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/out.pov",
                    required=False, dest="output",
                    help="output filename")

args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
headerFile = args.head
output_filename = args.output

max_z = 350
shift = 10

# Read the trj and gro file
u = Universe(psf_filename, traj_filename)

# Obtain initial information form gro and trj files
print "Total number of atoms = " + str(len(u.atoms))
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames)
#frame = num_frames / 2
frame = 1
f1 = open('/home/ramezani/midway/lc/' + headerFile, 'r')
f2 = open(output_filename, 'w+')
head = f1.read()
f1.close()
f2.write(head)

for curr_frame in xrange(0, frame + 1):
    if curr_frame != 0:
        trj = u.trajectory.next()
    else:
        trj = u.trajectory[0]
    if curr_frame < frame:
        continue
    else:
        box = u.dimensions[0:3]
        coor = u.selectAtoms("resname 5CB").coordinates()
        coor_N = u.selectAtoms("resname 5CB and name NY1").coordinates()
        coor_C = u.selectAtoms("resname 5CB and name CA12").coordinates()
        f2.write('// 5CB \n')
        R = 0.5
        G = 0.5
        B = 0.5
        T = 1.0
        Ra = 1.0
        for i in xrange(0, len(coor) - 19, 19):
            flag = False
            for j in range(0, 19):
                pos1 = coor[i + j]
                if pos1[2] > max_z:
                    flag = True
                    break
            if flag:
                for j in range(0, 19):
                    coor[i + j][2] = coor[i+j][2] - box[2]
                flag = False

        for i in xrange(0,len(coor_N)) :
            pos1 = coor[i*19]
            f2.write('sphere { < %5.3f, %5.3f, %5.3f>, %5.3f pigment {color rgb<%5.3f, %5.3f, %5.3f>} finish {phong .8 }} \n' % (pos1[0],pos1[1],pos1[2],Ra,0,0,1.))
            for j in range(i*19+1,i*19+19) :
                pos1 = coor[j]
                f2.write('sphere { < %5.3f, %5.3f, %5.3f>, %5.3f pigment {color rgb<%5.3f, %5.3f, %5.3f>} finish {phong .8 }} \n' % (pos1[0],pos1[1],pos1[2],Ra,R,G,B))

        # Write the pov-ray file for F9
        coor_f = u.selectAtoms("resname F9").coordinates()
        f2.write('// F \n')

        R = 1.0
        G = 0.6
        B = 0.6
        r = 1.0

        for i in xrange(0, len(coor_f) - 29, 29):
            flag = False
            for j in range(0, 29):
                pos1 = coor_f[i + j]
                if pos1[2] > max_z:
                    flag = True
                    break
            if flag:
                for j in range(0, 29):
                    coor_f[i + j][2] = coor_f[i+j][2] - box[2]
                flag = False

        for i in xrange(0, len(coor_f)):
            pos1 = coor_f[i]
            f2.write(
                'sphere { < %5.3f, %5.3f, %5.3f>, %5.3f pigment {color rgb<%5.3f, %5.3f, %5.3f>} finish {specular 0.7 roughness 0.03}} \n' % (
                pos1[0], pos1[1], pos1[2], r, R, G, B))

f2.close()

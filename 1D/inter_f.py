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
parser.add_argument("-l",
                    action="store", nargs='?', default=float(0),
                    required=False, dest="lowerBond", type=float,
                    help="the lower bond of your scale bar")
parser.add_argument("-u",
                    action="store", nargs='?', default=float(1),
                    required=False, dest="upperBond", type=float,
                    help="the upper bond of your scale bar")
parser.add_argument("-r",
                    action="store", nargs='?', default="head",
                    required=False, dest="head",
                    help="the upper bond of your scale bar")
parser.add_argument("-o",
                    action="store", nargs='?', default="../Analysis/out.pov",
                    required=False, dest="output",
                    help="output filename")
parser.add_argument("-i",
                    action="store", nargs='?', default="../Analysis/ScalerZ.dat",
                    required=False, dest="input",
                    help="input` filename")

args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
colorLow = args.lowerBond
colorHigh = args.upperBond
headerFile = args.head
output_filename = args.output
input_filename = args.input
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

f1 = open(input_filename, 'r')
z = []
S = []
for line in f1:
    line = line.strip()
    if line != "":
        columns = line.split()
        if '#' in columns[0]:
            check = 'yes'
        else:
            zz = float(columns[0])
            SS = float(columns[1])
            S.append(SS)
            z.append(zz)

y = np.zeros((len(z)), dtype='float')
for ii in range(0, len(z)):
    y[ii] = S[ii]

zmin = z[0]
zmax = z[len(z) - 1]
delta = 2
f = interpolate.interp1d(z, y, kind='linear')
znew = np.arange(zmin, zmax, delta)
ynew = f(znew)
# plt.plot(z, y, 'o', znew, ynew, '-')
# plt.show()

R = 0.0
G = 0.0
B = 0.0
T = 1.0
Ra = 1.0
for curr_frame in xrange(0, frame + 1):
    if curr_frame != 0:
        trj = u.trajectory.next()
    else:
        trj = u.trajectory[0]
    if curr_frame < frame:
        continue
    else:
        box = u.dimensions[0:3]
        coor = u.selectAtoms("name NY1").residues.coordinates()

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

        for i in xrange(0, len(coor)):
            pos1 = coor[i]
            nz = int((pos1[2] - zmin - delta / 2) / delta)
            if nz >= len(ynew):
                nz = len(ynew) - 1
            if nz < 0:
                nz = 0
            yy = ynew[nz]
            yy = (yy - colorLow) / (colorHigh - colorLow)

            B = exp(-(yy * yy) / 0.1)
            G = exp(-(yy - 0.5) * (yy - 0.5) / 0.1)
            R = exp(-(yy - 1.0) * (yy - 1.0) / 0.1)

            VV = sqrt(R * R + G * G + B * B)
            R /= VV
            G /= VV
            B /= VV
            R += 0.5
            G += 0.5
            B += 0.5

            f2.write(
                'sphere { < %5.3f, %5.3f, %5.3f>, %5.3f pigment {color rgb<%5.3f, %5.3f, %5.3f>} finish {phong .8 }} \n' % (
                pos1[0], pos1[1], pos1[2], Ra, R, G, B))
        # Write the pov-ray file for F9
        coor_f_c = u.selectAtoms("type CCF3 or type CCF2").coordinates()
        coor_f_f = u.selectAtoms("type FCF3 or type FCF2").coordinates()
        coor_f = u.selectAtoms("resname F9").coordinates()
        f2.write('// F \n')

        R = 1.0
        G = 0.6
        B = 0.6
        r = 1.0

        for i in xrange(0, len(coor_f) - 29 +1, 29):
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

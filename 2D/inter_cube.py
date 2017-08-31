#! /usr/bin/env python

from MDAnalysis import *
from math import *
import argparse
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Prepare POV-Ray input file for rendering")

parser.add_argument("-t",
                    action="store", nargs='?',
                    required=True, dest="traj",
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-l",
                    action="store", nargs='?', default= float(0),
                    required=False, dest="lowerBond", type = float,
                    help="the lower bond of your scale bar")
parser.add_argument("-u",
                    action="store", nargs='?', default = float(1),
                    required=False, dest="upperBond", type = float,
                    help="the upper bond of your scale bar")
parser.add_argument("-r",
                    action="store", nargs='?', default = "head",
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
parser.add_argument("-c",
                    action="store", nargs='?', default = int(2),
                    required=False, dest="columnNum", type = int,
                    help="the column number in the output file for plotting")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
colorLow = args.lowerBond
colorHigh = args.upperBond
headerFile = args.head
output_filename = args.output
input_filename = args.input
columnNum = args.columnNum

print traj_filename
# Read the trj and gro file
u = Universe(psf_filename, traj_filename)
box_v = u.dimensions[0:3]
box   = box_v[0]/2
# Obtain initial information form gro and trj files
print "Total number of atoms = " + str(len(u.atoms))
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames)
#frame = num_frames/2
frame = 1

f1=open('/home/ramezani/midway/lc/CUBE/'+headerFile,'r')
f2=open(output_filename, 'w+')
head = f1.read()
f1.close()
f2.write(head)

f1=open(input_filename,'r')
binsX=[]
binsY=[]
scalar=[]
mesh = []
x=0
for line in f1:
    line = line.strip()
    if line != "":
        columns = line.split()
        if '#' in columns[0]:
            check='yes'
        else:
            firstCol = float(columns[0])
            binsX.append(firstCol)
            mesh.append(firstCol)
            secCol = float(columns[1])
            binsY.append(secCol)
            mesh.append(secCol)
            thirdCol = float(columns[columnNum])
            scalar.append(thirdCol)

    x +=1

print 'Finished with reading data'

binsX=np.array(binsX)
binsY=np.array(binsY)
scalar=np.array(scalar)
scalar=np.reshape(scalar,(np.sqrt(len(binsX)),np.sqrt(len(binsY))))

nn = np.sqrt(len(binsX))
x = np.mgrid[-box:box:complex(nn)]
y = np.mgrid[-box:box:complex(nn)]
f = interpolate.interp2d(x, y, scalar, kind='cubic')

delta =2
xnew = np.arange(-box, box, delta)
ynew = np.arange(-box, box, delta)

znew = f(xnew, ynew)


R = 0.0
G = 0.0
B = 0.0
T = 1.0
Ra = 1.0
for curr_frame in xrange(0, frame+1) :
    if curr_frame != 0 :
        trj = u.trajectory.next()
    else :
        trj = u.trajectory[0]
    if curr_frame < frame :
        continue
    else :
#        coor = u.selectAtoms("resname 5CB").coordinates()
        coor = u.selectAtoms("all").coordinates()
        for i in xrange(0,len(coor)) :
            pos1 = coor[i]
            nx = int((pos1[0]+box)/delta)
            if nx>=len(xnew):
                nx=len(xnew)-1
            if nx<0 :
                nx=0
            ny = int((pos1[1]+box)/delta)
            if ny>=len(ynew):
                ny=len(ynew)-1
            if ny<0 :
                ny=0

            zz = znew[nx,ny]
            zz = (zz-colorLow)/(colorHigh - colorLow)

            R = exp(-(zz*zz)/0.1)
            G = exp(-(zz-0.5)*(zz-0.5)/0.1)
            B = exp(-(zz-1.0)*(zz-1.0)/0.1)

            VV = sqrt(R*R+G*G+B*B)

            if VV != 0:
                R /= VV
                G /= VV
                B /= VV
                R += 0.5
                G += 0.5
                B += 0.5

                f2.write('sphere { < %5.3f, %5.3f, %5.3f>, %5.3f pigment {color rgb<%5.3f, %5.3f, %5.3f>} finish {phong .8 }} \n' % (pos1[0],pos1[1],pos1[2],Ra,R,G,B))

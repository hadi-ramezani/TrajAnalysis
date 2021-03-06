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
from scipy import interpolate
from scipy.optimize import curve_fit
import scipy.integrate


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
                    help="Number of atoms in each slab")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")
parser.add_argument("-o1",
                    action="store", nargs='?', default="../Analysis/anchoring.dat",
                    required=False, dest="output1",
                    help="output filename for Anchoring")
parser.add_argument("-o2",
                    action="store", nargs='?', default="../Analysis/hist_P.dat",
                    required=False, dest="output2",
                    help="output filename for Hist_P")
parser.add_argument("-o3",
                    action="store", nargs='?', default="../Analysis/hist_S.dat",
                    required=False, dest="output3",
                    help="output filename for Hist_S")
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
parser.add_argument("-T",
                    action="store", nargs='?',
                    required=True, dest="temperature", type=float,
                    help="Simulation temperature")
parser.add_argument("-sn",
                    action="store", nargs='?', default="name NY1",
                    required=False, dest="atom_selection_N",
                    help="Atom selection argument for N")
parser.add_argument("-sc",
                    action="store", nargs='?', default="name CA12",
                    required=False, dest="atom_selection_C",
                    help="Atom selection argument for C")
parser.add_argument("-s",
                    action="store", nargs='?', default="circle",
                    required=False, dest="shape",
                    help="shape of the bubble: circle, cube, squircle")



args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
nslab = args.nslab
output_filename1 = args.output1
output_filename2 = args.output2
output_filename3 = args.output3
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
temperature = args.temperature
atom_selection_N = args.atom_selection_N
atom_selection_C = args.atom_selection_C
shape = args.shape
time1 = time.time()
output_file1 = open(output_filename1,"w")
output_file2 = open(output_filename2,"w")
output_file3 = open(output_filename3,"w")


def distance(v1, v2) :
    return sqrt((v1[0] - v2[0])*(v1[0] - v2[0]) + (v1[1] - v2[1])*(v1[1] - v2[1]) + (v1[2] - v2[2])*(v1[2] - v2[2]) )
def distance_d(v1, v2) :
    return abs(v1-v2)
def vector(v1, v2):
    return [v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]
def position(v1,v2):
    posx = (v1[0]+v2[0])/2.0
    posy = (v1[1]+v2[1])/2.0
    if shape == "circle":
        return np.sqrt(posx*posx+posy*posy)
    elif shape == "cube":
        return np.power(np.power(posx,40)+np.power(posy,40),0.025)
    elif shape == "squircle":
        return np.power(np.power(posx,4)+np.power(posy,4),0.25)
    else:
        print sys.stderr,"Unknown shape, exiting ..."
        sys.exit(1)

def Sz(v) :
    return (3.*v[2]*v[2]-1.)/2.0
def cosZ(v) :
	return v[2]
def get_director(loc_C, loc_N):
    ''' This function calculate the director using gradient of the
     shape'''
    n = np.zeros(3)
    if shape == "circle":
        n[0] = (loc_C[0]+loc_N[0])/2.0
        n[1] = (loc_C[1]+loc_N[1])/2.0
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    elif shape == "cube":
        n[0] = np.power((loc_C[0]+loc_N[0])/2.0,39)
        n[1] = np.power((loc_C[1]+loc_N[1])/2.0,39)
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    elif shape == "squircle":
        n[0] = np.power((loc_C[0]+loc_N[0])/2.0,3)
        n[1] = np.power((loc_C[1]+loc_N[1])/2.0,3)
        n[2] = 0.0
        n = n/np.linalg.norm(n)
    else:
        print sys.stderr,"Unknown shape, exiting ..."
        sys.exit(1)
    return n


def get_area(binCenter):
    if shape == "circle":
        A = 2*pi*binCenter*box[2]
        #A = np.zeros(len(binCenter))
        #for i in range(len(binCenter)):
        #    A[i] , B = scipy.integrate.quad(lambda x : np.sqrt(1+np.power(x,2)/(np.power(binCenter[i],2)-np.power(x,2))), -binCenter[i], binCenter[i])
        #    A[i] = 2 * A[i] * box[2]
    elif shape == "cube":
        A = np.zeros(len(binCenter))
        for i in range(len(binCenter)):
            # Calculating the perimeter of cube using integral(-R,R)(sqrt(1+(dy/dx)^2)dx (numerical integration). Only y>0 is considered first. The calculated value is
            # the multiplied by 2 to yield perimeter of the shape. Multiplication of the perimeter by z dimension of the box gives the total area.
            # For more detail see http://www.mathwords.com/a/arc_length_of_a_curve.htm
            A[i] , B = scipy.integrate.quad(lambda x : np.sqrt(1+np.power(x,78)/np.power((np.power(binCenter[i],40)-np.power(x,40)),1.95)), -binCenter[i], binCenter[i])
            A[i] = 2 * A[i] * box[2]
    elif shape == "squircle":
        A = np.zeros(len(binCenter))
        for i in range(len(binCenter)):
            A[i] , B = scipy.integrate.quad(lambda x : np.sqrt(1+np.power(x,6)/np.power((np.power(binCenter[i],4)-np.power(x,4)),1.5)), -binCenter[i], binCenter[i])
            A[i] = 2 * A[i] * box[2]
    return A


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
#delta = nslab
#delta = 50
print "Number of atoms in each slab is around " + str(delta)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames)
if end_frame > num_frames :
  end_frame = num_frames
print "reading frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays

bins = np.zeros((nslab+1),dtype='float')

pvector = np.zeros((natoms),dtype='float')

hist_cos_P = np.empty((0),dtype='float')
hist_cos_S = np.empty((0),dtype='float')

dens = np.zeros((nslab),dtype='int')
dens_f = np.empty((0),dtype='float')

# Boltzmann constant and kT
k = 0.0013806488
kT = -k * temperature


nframes = 0
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

        # Applay PBC 
        box = u.dimensions[0:3]

        #Calculate the position of each slab to have almost same number of atoms in each slab
#        if curr_frame==start_frame :
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            pvector[i] = position(loc_C[0:2],loc_N[0:2])

        pvector = np.sort(pvector)
        bins[nslab] = pvector[natoms-1]
        for i in range(nslab):
            bins[i] = pvector[i*delta]
#        print bins

        dens.fill(0)
        hist_cos_P_temp = np.empty((0))
        hist_cos_S_temp = np.empty((0))

        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            p = position(loc_C[0:2],loc_N[0:2])
            x = 1
            while x<=nslab :
              if p<=bins[x] :
                 vv = vector(loc_C,loc_N)
                 vv = vv/np.linalg.norm(vv)
                 n  = get_director(loc_C,loc_N)
                 cos_S = np.dot(vv,n)
                 hist_cos_S_temp = np.append(hist_cos_S_temp,p)
                 hist_cos_S_temp = np.append(hist_cos_S_temp,cos_S)
                 dens[x-1] += 1

                 break
              else :
                 x += 1
        hist_cos_S = np.append(hist_cos_S,hist_cos_S_temp,0)
        dens_f     = np.append(dens_f,dens,0)

# Calculate average density
dens_f       = np.reshape(dens_f,( len(dens_f)/nslab,nslab))
dens_average =  np.average(dens_f,0)

hist_cos_S = np.reshape(hist_cos_S, (len(hist_cos_S)/2,2))
hist_cos_P = hist_cos_S

# Define bins of cos(beta)
cosGrid = np.linspace(-1,1,num=21)

# Histogram all data
H, yedges, xedges = np.histogram2d(hist_cos_P[:,0], hist_cos_P[:,1], bins=(bins,cosGrid))

P_H = np.empty(H.shape)

# Define bin centers
binCenter     = np.empty(len(bins)-1)
cosGridCenter = np.empty(len(cosGrid)-1)

# Calculate bin center
for x in range(len(cosGridCenter)):
    cosGridCenter[x] = (cosGrid[x]+cosGrid[x+1])/2

# Calculate probability from histogram data
for x in range(nslab):
    binCenter[x]     = (bins[x]+bins[x+1])/2
    sum = np.sum(H[x,:])
    P_H[x,:] = H[x,:]/sum

# Write out the hist_P data
for x in range(nslab):
    for y in range(len(cosGridCenter)):
        output_file2.write('%5.3f %7.5f %7.5f\n' % ( binCenter[x] ,cosGridCenter[y], P_H[x,y]))

# Calculate surface free energy and calculate input data for anchoring strength
F = np.empty(P_H.shape)
xdata = np.empty(P_H.shape)
sigma = np.empty(P_H.shape)

A  = get_area(binCenter)

x0 = np.array([0.001,0.001])

for i in range(nslab):
    for j in range(len(cosGridCenter)):
        if P_H[i,j] == 0 :
            P_H[i,j] = 1e-8

        sigma[i,j] = 1/sqrt(P_H[i,j])
        F[i,j]= kT * log(P_H[i,j]) * dens[i]/A[i]

        beta = acos(cosGridCenter[j])
        betaEq = 0
        xdata[i,j] = pow(cos(beta - betaEq),2)

# Calculate anchoring strength
def func(x, W0, W2):
    return W0 - 0.5 * W2 * x

popt = np.zeros((nslab,2))

# Do the optimization to obtain W0 and W2
for i in range(nslab):
    popt[i,:], pcov = curve_fit(func, xdata[i,:], F[i,:], p0=x0, sigma=sigma[i,:])

# Write out the anchoring strength data
for x in range(nslab):
    output_file1.write('%5.3f %7.5f %7.5f\n' % ( binCenter[x] ,popt[x,0], popt[x,1]))

# Plot the anchoring data
x = 0

# Get the limit of the plot based on distance from the surface
for x in range(nslab):
    if (binCenter[x] - binCenter[0]) < 50:
        x1=x

    elif (binCenter[nslab-1] - binCenter[x]) < 50:
        x2=x
        break

plt.plot(binCenter[0:x1], popt[0:x1,1], marker='o', color='b', markeredgecolor='b')
#plt.plot(binCenter[x2:], popt[x2:,1], marker='o', color='b', markeredgecolor='b')

plt.rcParams['axes.linewidth'] = 1.5
plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.ylabel(r'$\mathrm{W_2^A(J/m^2)}$', fontsize=27)
plt.tick_params(which='both', width=2)
plt.xlim((binCenter[0],binCenter[x1]))
#plt.ylim((0,0.05))
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(10,5)
plt.savefig('../Analysis/anchoring.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

plt.clf()

# Do the interpolation for distribution graphs
xi = np.mgrid[-1:1:500j]
yi = np.mgrid[binCenter[0]:binCenter[len(binCenter)-1]:1000j]
f = interpolate.interp2d(cosGridCenter, binCenter, P_H, kind='cubic')
zi = f(xi, yi)

# Plot hist_cos_P data
plt.gcf().set_size_inches(7,10)
plt.xlim(-1,1)
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))
plt.contour(xi,yi,zi, linewidths=1, colors='k')
plt.pcolormesh(xi,yi,zi,  vmin=0, vmax=0.2, cmap='Reds' )
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$\mathrm{P(z,cos\beta)}$', fontsize=17)
plt.xlabel(r'$\mathrm{cos\beta}$', fontsize=27)
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.tick_params(which='both', width=2)
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.rcParams['axes.linewidth'] = 1.5
plt.savefig('../Analysis/plot_P_cos_his.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
plt.clf()

# Plot hist_cos_S data
cosGrid = np.linspace(-1,1,num=20)
H, yedges, xedges = np.histogram2d(hist_cos_S[:,0], hist_cos_S[:,1], bins=(bins,cosGrid))
P_H = np.empty(H.shape)
binCenter = np.empty(len(bins)-1)
cosGridCenter = np.empty(len(cosGrid)-1)
for x in range(len(cosGridCenter)):
    cosGridCenter[x] = (cosGrid[x]+cosGrid[x+1])/2

for x in range(nslab):
    binCenter[x]     = (bins[x]+bins[x+1])/2
    sum = np.sum(H[x,:])
    P_H[x,:] = H[x,:]/sum

# Write out the hist_S data
for x in range(nslab):
    for y in range(len(cosGridCenter)):
        output_file3.write('%5.3f %7.5f %7.5f\n' % ( binCenter[x] ,cosGridCenter[y], P_H[x,y]))

xi = np.mgrid[-1:1:500j]
yi = np.mgrid[binCenter[0]:binCenter[len(binCenter)-1]:1000j]
f = interpolate.interp2d(cosGridCenter, binCenter, P_H, kind='cubic')
zi = f(xi, yi)
plt.gcf().set_size_inches(7,10)
plt.xlim(-1,1)
plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))
plt.contour(xi,yi,zi, linewidths=1, colors='k')
plt.pcolormesh(xi,yi,zi,  vmin=0, vmax=0.2, cmap='Reds' )
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$\mathrm{P(z,cos\beta)}$', fontsize=17)
plt.xlabel(r'$\mathrm{cos\beta}$', fontsize=27)
plt.ylabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.tick_params(which='both', width=2)
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.rcParams['axes.linewidth'] = 1.5
plt.savefig('../Analysis/plot_S_cos_his.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
plt.clf()

time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file1.close()
output_file2.close()
output_file3.close()
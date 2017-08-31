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
from scipy.interpolate import spline


#Plot the output results
outFile = open('../Analysis/dist_N.dat', 'r')
z = []
P = []
for line in outFile:
	line = line.strip()
	if line != "":
		columns = line.split()
		if '#' in columns[0]:
			check='yes'
		else:
			firstCol = float(columns[0])
			z.append(firstCol)
			secCol = float(columns[1])
			P.append(secCol) 

#Plot the P2 data
z_smooth = np.linspace(z[0], z[len(z)-1], len(z))
P_smooth = spline(z, P, z_smooth, 3)
plt.plot(z_smooth, P_smooth, linewidth=2, color='b' , label="LC N")
plt.hold(True)
###################################################################################################
###################################################################################################
outFile = open('../Analysis/dist_O.dat', 'r')
z = []
P = []
for line in outFile:
	line = line.strip()
	if line != "":
		columns = line.split()
		if '#' in columns[0]:
			check='yes'
		else:
			firstCol = float(columns[0])
			z.append(firstCol)
			secCol = float(columns[1])
			P.append(secCol) 

#Plot the P2 data
z_smooth = np.linspace(z[0], z[len(z)-1], len(z))
P_smooth = spline(z, P, z_smooth, 3)
plt.plot(z_smooth, P_smooth, linewidth=2, color='r' , label="Water O")
plt.hold(True)
###################################################################################################
###################################################################################################
outFile = open('../Analysis/dist_cation.dat', 'r')
z = []
P = []
for line in outFile:
	line = line.strip()
	if line != "":
		columns = line.split()
		if '#' in columns[0]:
			check='yes'
		else:
			firstCol = float(columns[0])
			z.append(firstCol)
			secCol = float(columns[1])
			P.append(secCol) 

#Plot the P2 data
z_smooth = np.linspace(z[0], z[len(z)-1], len(z))
P_smooth = spline(z, P, z_smooth, 3)
plt.plot(z_smooth, P_smooth, linewidth=2, color='y' , label="Na")
plt.hold(True)
###################################################################################################
###################################################################################################
outFile = open('../Analysis/dist_anion.dat', 'r')
z = []
P = []
for line in outFile:
	line = line.strip()
	if line != "":
		columns = line.split()
		if '#' in columns[0]:
			check='yes'
		else:
			firstCol = float(columns[0])
			z.append(firstCol)
			secCol = float(columns[1])
			P.append(secCol) 

#Plot the P2 data
z_smooth = np.linspace(z[0], z[len(z)-1], len(z))
P_smooth = spline(z, P, z_smooth, 3)
plt.plot(z_smooth, P_smooth, linewidth=2, color='g' , label="Cl")
plt.hold(True)
###################################################################################################
###################################################################################################
# Legend properties
plt.rcParams['axes.linewidth'] = 1.5
plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)
plt.ylabel(r'$\mathrm{P(z)}$', fontsize=27)	
plt.tick_params(which='both', width=2)		
plt.xlim((0,350))
plt.ylim((0,0.03))       
plt.yticks(fontsize=17)
plt.xticks(fontsize=17)
plt.gcf().set_size_inches(10,5)
plt.legend(loc='best', frameon=True, fancybox=True)
plt.savefig('../Analysis/dist.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
	

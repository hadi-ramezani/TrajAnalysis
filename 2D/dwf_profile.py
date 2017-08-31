#! /usr/bin/env python

from MDAnalysis import *
from math import *
import argparse
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot density profile of cube-like systems"
                                             " in different directions")

parser.add_argument("-t",
                    action="store", nargs='?',
                    required=True, dest="traj",
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-o",
                    action="store", nargs='?', default="density_profile",
                    required=False, dest="output",
                    help="output filename")
parser.add_argument("-i",
                    action="store", nargs='?', default="../Analysis/dwf.dat",
                    required=False, dest="input",
                    help="input filename")
parser.add_argument("-s",
                    action="store", nargs='?', default="circle",
                    required=False, dest="shape",
                    help="shape of the bubble: circle, cube, squircle")
parser.add_argument("-c",
                    action="store", nargs='?', default = int(2),
                    required=False, dest="columnNum", type = int,
                    help="the column number in the output file for plotting")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
output_filename = args.output
input_filename = args.input
shape = args.shape
columnNum = args.columnNum

print traj_filename
# Read the trj and dcd file
u = Universe(psf_filename, traj_filename)
box_v = u.dimensions[0:3]
box   = box_v[0]/2
def read_S(input_filename):
    # Read S data from the input file
    inputFile = open(input_filename,'r')
    binsX = []
    binsY = []
    scalar = []
    mesh = []
    density = []
    x = 0
    for line in inputFile:
        line = line.strip()
        if line != "":
            columns = line.split()
            if not line.startswith("#"):
                binsX.append(float(columns[0]))
                mesh.append(float(columns[0]))
                binsY.append(float(columns[1]))
                mesh.append(float(columns[1]))
                scalar.append(float(columns[columnNum]))
                density.append(float(columns[4]))
        x += 1

    print 'Finished with reading data'
    binsX = np.array(binsX)
    binsY = np.array(binsY)
    mesh = np.array(mesh)
    scalar = np.array(scalar)
    density = np.array(density)
    scalar_reshaped = np.reshape(scalar,(np.sqrt(len(binsX)),np.sqrt(len(binsY))))
    density_reshaped = np.reshape(density,(np.sqrt(len(binsX)),np.sqrt(len(binsY))))
    mesh = np.reshape(mesh,(len(mesh)/2,2))
    return binsX, binsY, scalar, scalar_reshaped, density_reshaped


def intepolate_S(binsX, scalar_reshaped, delta):
    # Do the interpolation to increase the resolution
    nn = np.sqrt(len(binsX))
    x = np.mgrid[-box:box:complex(nn)]
    y = np.mgrid[-box:box:complex(nn)]
    f = interpolate.interp2d(x, y, scalar_reshaped, kind='cubic')

    delta = 1
    num = 2 * box /delta
    xnew = np.linspace(-box, box, num)
    ynew = np.linspace(-box, box, num)
    znew = f(xnew, ynew)
    return xnew, ynew, znew

def get_r(x, y):
    if shape == "circle":
        return np.sqrt(x*x+y*y)
    elif shape == "cube":
        return np.power(np.power(x,40)+np.power(y,40),0.025)
    elif shape == "squircle":
        return np.power(np.power(x,4)+np.power(y,4),0.25)
    else:
        print sys.stderr,"Unknown shape, exiting ..."
        sys.exit(1)

def s_density_axis(xnew, ynew, znew):
    length = len(xnew)
    # Find the S values along secondary diagonal
    x_sec_dia = []
    y_sec_dia = []
    scalar_sec_dia = []
    for i, x in enumerate(xnew):
        for j, y in enumerate(ynew):
            if int(x) == int(y) and not abs(int(x)) in x_sec_dia:
                x_sec_dia.append(abs(x))
                y_sec_dia.append(abs(y))
                scalar_sec_dia.append((znew[i,j] + znew[length-i-1,length-j-1]) / 2)
    x_sec_dia = np.array(x_sec_dia); y_sec_dia = np.array(y_sec_dia)
    scalar_sec_dia = np.array(scalar_sec_dia)
    r_sec_dia = get_r(x_sec_dia,y_sec_dia)

    # Find the S values along primary diagonal
    x_prim_dia = []
    y_prim_dia = []
    scalar_prim_dia = []
    for i, x in enumerate(xnew):
        for j, y in enumerate(ynew):
            if int(x) == int(-y) and not abs(int(x)) in x_prim_dia:
                x_prim_dia.append(abs(x))
                y_prim_dia.append(abs(y))
                scalar_prim_dia.append((znew[i,j] + znew[length-i-1,length-j-1] + znew[j,i] + znew[length-j-1,length-i-1]) / 4)
    x_prim_dia = np.array(x_prim_dia); y_prim_dia = np.array(y_prim_dia)
    scalar_prim_dia = np.array(scalar_prim_dia)
    r_prim_dia = get_r(x_prim_dia, y_prim_dia)

    # Find the S values along y = 0 or x = 0
    x_xy_0 = []
    y_xy_0 = []
    scalar_x_0 = []
    scalar_y_0 = []
    for i, x in enumerate(xnew):
        for j, y in enumerate(ynew):
            if int(x) == 0 and x >= 0:
                x_xy_0.append(abs(x))
                y_xy_0.append(abs(y))
                scalar_x_0.append((znew[i,j]+znew[length-i-1, length-j-1]) / 2)
            if int(y) == 0 and y >= 0:
                scalar_y_0.append((znew[i,j]+znew[length-i-1, length-j-1]) / 2)
    x_xy_0 = np.array(x_xy_0); y_xy_0 = np.array(y_xy_0)
    scalar_x_0 = np.array(scalar_x_0); scalar_y_0 = np.array(scalar_y_0)
    scalar_xy_0 = (scalar_x_0 + scalar_y_0) / 2
    r_xy_0 = get_r(x_xy_0, y_xy_0)

    return r_sec_dia, scalar_sec_dia, r_prim_dia, scalar_prim_dia, r_xy_0, scalar_xy_0

def plot_data(r_sec_dia, s_dens_sec_dia, r_prim_dia, s_dens_prim_dia, r_xy_0,
              s_dens_xy_0, outputname, ylabel_txt):
    # create the general figure
    fig1 = plt.figure()
    ax1 = plt.subplot(111)
    plt.subplots_adjust(hspace=0.001)
    plt.gcf().set_size_inches(6,5)

    ax1.plot(r_prim_dia, s_dens_prim_dia, marker='o', color='cyan' , label="Main Dia.")
    ax1.plot(r_sec_dia, s_dens_sec_dia, marker='o', color='fuchsia' , label="Sec. Dia.")
    ax1.plot(r_xy_0, s_dens_xy_0, marker='o', color='orange' , label="x=0 or y=0")

    plt.legend(loc='best', frameon=True, fancybox=True, prop={'size':10})
    plt.xlabel(r'$\mathrm{r\/(\AA)}$', fontsize=25, labelpad=0)
    plt.ylabel(ylabel_txt, fontsize=20)
    plt.tick_params(which='both', width=2)
    plt.xlim(50,160)
    plt.ylim((0, 15))
    #ax1.yaxis.set_ticks(np.arange(-18, 3, 3))
    # plt.yticks(weight='bold')
    # plt.xticks(weight='bold')
    plt.yticks(fontsize=17)
    plt.xticks(fontsize=17)
    plt.xticks(rotation=20)
    plt.savefig(outputname, dpi=1000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)

# Calculate and plot scalar order profiles
binsX, binsY, scalar, scalar_reshaped, density_reshaped = read_S(input_filename)
xnew, ynew, znew = intepolate_S(binsX, scalar_reshaped, delta=1)
r_sec_dia, scalar_sec_dia, r_prim_dia, scalar_prim_dia, r_xy_0, scalar_xy_0 = s_density_axis(xnew, ynew, znew)

plot_data(r_sec_dia[:len(r_sec_dia)/2], scalar_sec_dia[:len(scalar_sec_dia)/2],
          r_prim_dia[:len(r_prim_dia)/2], scalar_prim_dia[:len(scalar_prim_dia)/2],
          r_xy_0[:len(r_xy_0)/2], scalar_xy_0[:len(scalar_xy_0)/2],
          outputname = '../Analysis/'+ output_filename+'.png',
          ylabel_txt = r'$\mathrm{\langle u^{2}\rangle\/(\AA^{2})}$')

# Calculate and plot number density profiles
'''xnew_density, ynew_density, znew_density = intepolate_S(binsX, density_reshaped, delta=1)
r_sec_dia_dens, density_sec_dia, r_prim_dia_dens, density_prim_dia, r_xy_0_dens, density_xy_0 = s_density_axis(xnew, ynew, znew_density)

plot_data(r_sec_dia_dens[:len(r_sec_dia_dens)/2], density_sec_dia[:len(density_sec_dia)/2],
          r_prim_dia_dens[:len(r_prim_dia_dens)/2], density_prim_dia[:len(density_prim_dia)/2],
          r_xy_0_dens[:len(r_xy_0_dens)/2], density_xy_0[:len(density_xy_0)/2],
          outputname = '../Analysis/density_axis.png',
          ylabel_txt = r'$\mathrm{{density\/(\frac{mol.}{nm^{3}})}}$')'''

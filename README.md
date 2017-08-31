# TrajAnalysis

This project contains several python scripts that I wrote over years to analyze the 
trajectories of atomistic MD simulations. The scripts make use 
of MDAnalysis library (http://www.mdanalysis.org/). I used version 0.10.0 of the 
library at the time. The syntax might be slightly different in newer versions.
Most of the scripts should be easily adaptable for use in similar systems such 
as lipid membranes.  Some scripts perform trivial calculations and some others
accomplish pretty sophisticated analysis. Currently, there are two main folders
in this project.

1) 1D: This folder contains scripts to analyze thin films of liquid crystals 
(free-standing films and systems in contact with aqueous phases). Please 
cite the following papers when you use the scripts in this folder:

- “Understanding atomic-scale behavior of liquid crystals at aqueous interfaces”
H. Ramezani-Dakhel, M. Sadati, M. Rahimi, A. Ramirez-Hernandez, Benoît Roux, and
J. J. de Pablo Journal of Chemical Theory and Computation, 13, 237-244, 2017.

- “Molecular structure of canonical liquid crystal interfaces”
M. Sadati,* H. Ramezani-Dakhel,* W. Bu, E. Sevgen, Z. Liang, C. Erol, 
N. Taheri Qazvini, M. Rahimi, B. Lin, N. L. Abbott, B. Roux, M. Schlossman,
J. J. de Pablo, Journal of the American Chemical Society, 139, 3841-3850, 2017.

2) 2D: This folder contains scripts for two-dimensional analysis of liquid 
crystalline systems. The systems of interests will be non-uniform in x and y 
directions and uniform in z direction. Please cite the following paper when you 
use these scripts:

- “Segregation of liquid crystals mixture in topological defects”
M. Rahimi, H. Ramezani-Dakhel, R. Zhang, A. Ramirez-Hernandez, N. L. Abbott, 
and J. J. de Pablo, Nature Communications, 8, 15064, 2017.

Here I provide a brief description of the main codes in this project. 
Please feel free to explore other files in the folders as you wish. I assume that you have
psf and dcd files from you simulations. In case you use other file formats, you 
still should be able to use the scripts with some small modifications.

1) 1D

Please note that all scripts can read a single trajectory file (e.g. myoutput.dcd)
or sequences of trajectory files that are named 1.dcd, 2.dcd, ..., N.dcd. In the
latter case you would need to provide the first and last numbers to the script, i.e.
1 and N. 

1.1. sZ_cosZ.py: This script computes the profile of P2 order parameter, number 
density, and average cos(beta), and plots them. An example command for running the 
script would be:

./sZ_cosZ.py  -p \<psf file> -t \<dcd file> -n \<number of bins> -b \<first frame> 
-e \<last frame>

Alternatively, you can do "./sZ_cosZ.py --help" to see all input options.

1.2. tensorZ.py: This script computes the profile of scalar order parameter (S), 
and components of nematic director (nk) and plots them for you.
To run the script you can do:

./tensorZ.py  -p \<psf file> -t \<dcd file> -n \<number of bins> -b \<first frame> 
-e \<last frame>

Alternatively, you can run "./tensorZ.py --help" to see all input options.

1.3. anchoring.py: This script calculates the profile of anchoring strength. To 
run the script you can issue the following command:

./anchoring.py  -p \<psf file> -t \<dcd file> -n \<number of bins> 
-T \<Temperature in Kelvin> -b \<first frame> -e \<last frame> 

Alternatively, you can do "./anchoring.py --help" to see all input options.

1.4. electron_density.py: This script computes the profile of electron density as 
a function of z. A sample command to run the script would be: 

./electron_density.py -p \<psf file> -t \<dcd file> -d \<bin size in angstrom> 
-s \<atom selection text> -o \<output file>  -b \<first frame> -e \<last frame>

Alternatively, you can run "./electron_density.py --help" to see all input options.

1.5. mass_density.py: Computes the mass density profile as a function of z. 
To run the script you can issue the following command:

./mass_density.py -p \<psf file> -t \<dcd file> -d \<bin size in angstrom> 
-s \<atom selection text> -o \<output file>  -b \<first frame> -e \<last frame>

Alternatively, you can just do "./mass_density.py --help" to see input options.

1.6. polarization.py: Calculates the polarization density profile of individual 
components in the system. The script currently works for a system that contains 
5CB, and water but it can easily be adapted to other systems. To use this script
you can do:

./polarization.py -p \<psf file> -t \<dcd file> -d \<bin size in angstrom> 
-b \<first frame> -e \<last frame>

Alternatively, you can just do "./polarization.py --help" to see all input options.

1.7. dielectric.py: This script computes dielectric constant profile of individual 
components in a thin film. To use the script you can do:

./dielectric.py  -p  \<psf filename> -t \<dcd filename> -d \<binsize in angstrom> 
-s \<atom selection command> -o \<output filename> -b \<first frame> -e \<last frame>

Alternatively, you can type "./dielectric.py --help" to see all input options.

1.8. dist.py: Computes distribution profile of individual components in a thin film.
An example command to run the script would be:

./dist.py -p \<psf filename> -t \<dcd filename> -d \<binsize in angstrom> 
-s \<atom selection command> -o \<output filename> -b \<first frame> -e \<last frame>

Alternatively, you can just do "./dist.py --help" to see all input options.

1.9. inter*.py: These scripts reads the output of other scripts (P2, S, etc), 
as well as the coordinates of atoms to generate inputs for POV-Ray. 
These are then used to make snapshots. 
An example command for using the script to generate the input for POV-Ray would be:

./inter*.py -t \<dcd filename> -p \<psf filename> -l \<lower color bound> 
-u \<upper color bound> -r \<header file> -o \<output name>   -i \<input name>

Alternatively, you can just do "./inter*.py --help" to see all input options.

2) 2D

Please note that all scripts can read a single trajectory file (e.g. myoutput.dcd)
or sequences of trajectory files that are named 1.dcd, 2.dcd, ..., N.dcd. In the
latter case you would need to provide the first and last number to the script, i.e.
1 and N. 

2.1. Q_vector.py: This script computes 2D map of scalar order parameter (S), 
and vectors of nematic director and outputs the data. To run the script you can do:

./Q_vector.py -p \<psf filename> -t \<dcd filename> 
-nx \<# of xbins for S calculations> -ny \<# of xbins for S calculations>
-vx \<# of xbins for director calculations> -vy \<# of xbins for director calculations>
-b \<first frame> -e \<last frame>

You can also do "Q_vector.py --help" to see more information.

2.2. densityProfile_defectArea.py: This script computes the profile of scalar
order parameter in different directions (e.g. main diagonal, secondary diagonal) 
as well as the size of the defect area in a two dimensional system. To use the 
script you can do:

./densityProfile_defectArea.py -p \<psf filename> -t \<dcd filename>
-s \<shape of the central bubble>

You can also do "densityProfile_defectArea.py --help" to see more information.

2.3. anchoring.py: This script calculates the profile of anchoring strength. To 
run the script you can issue the following command:

./anchoring.py  -p \<psf file> -t \<dcd file> -n \<number of bins> 
-T \<Temperature in Kelvin> -s \<shape of the central bubble> -b \<first frame> -e \<last frame>

You can also do "./anchoring.py --help" to see more information.

2.4. dwf.py: Compute a 2D map of Debye-Waller factor. To use the code you can do:

./dwf.py -p \<psf file> -t \<dcd file> -nx \<# of xbins> -ny \<# of xbins> 
-s \<shape of the central bubble> -b \<first frame> -e \<last frame>


You can also do "./dwf.py --help" to see more information.


********************************************************************************
Disclaimer- The scripts in this project are provided to you as is. Your feedbacks
are always welcome and I do my best to address them as soon as possible. However,
I accept no responsibilities about the accuracy of your calculations. It is your
own responsibility to ensure the accuracy and correctness of your results.
********************************************************************************

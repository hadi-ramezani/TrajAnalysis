#! /bin/bash

#lcDir="/project/roux/ramezani/lc/CUBE"
lcDir="/home/ramezani/midway/lc/CUBE"

#
###################################################################################################
#circle @ 25c
if true
then
	first=25
	last=25
	b=50
	e=100000
	#Navigate to the directory
	cd $lcDir/circle/DCD
	#Excecute the python script to calculate DWF
	#$lcDir/dwf.py -p ../w5cb_circle.psf -f $first -l $last -nx 100 -ny 100 -b $b -e $e -s circle
	#Execute the python script to prepare the POV-ray input files
	#$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_circle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_z.pov   -i ../Analysis/dwf.dat -c 4 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_circle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_sum.pov -i ../Analysis/dwf.dat -c 5 
	#$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_circle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_pr.pov   -i ../Analysis/dwf.dat -c 6 
	#$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_circle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_or.pov   -i ../Analysis/dwf.dat -c 7 
	# Excecute the python script to generate the dwf profiles along main directions
	#$lcDir/dwf_profile.py -p ../w5cb_circle.psf -t $last.dcd -s circle -o profile_dwf_z   -c 4
	$lcDir/dwf_profile.py -p ../w5cb_circle.psf -t $last.dcd -s circle -o profile_dwf_sum -c 5
	#$lcDir/dwf_profile.py -p ../w5cb_circle.psf -t $last.dcd -s circle -o profile_dwf_pr   -c 6
	#$lcDir/dwf_profile.py -p ../w5cb_circle.psf -t $last.dcd -s circle -o profile_dwf_or   -c 7
	#Navigate to the Analysis directory
	cd $lcDir/circle/Analysis
	#Plot the dwf data using gnuplot
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
	#Prepare the snapshots using POV-ray
	#povray +w3000 +H3000 -Iout_dwf_z.pov    +Osnap_dwf_z.png   +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_sum.pov  +Osnap_dwf_sum.png +Q9 +A0.01 +R9 -D -V
	#povray +w3000 +H3000 -Iout_dwf_pr.pov    +Osnap_dwf_pr.png   +Q9 +A0.01 +R9 -D -V
	#povray +w3000 +H3000 -Iout_dwf_or.pov    +Osnap_dwf_or.png   +Q9 +A0.01 +R9 -D -V
fi	
#
###################################################################################################
#cube @ 25c
if false
then
	first=3
	last=3
	b=5038
	e=100000
	#Navigate to the directory
	cd $lcDir/cube/DCD
	#Excecute the python script to calculate DWF
	#$lcDir/dwf.py -p ../w5cb_cube.psf -f $first -l $last -nx 100 -ny 100 -b $b -e $e -s cube
	#Execute the python script to prepare the POV-ray input files
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_cube.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_z.pov   -i ../Analysis/dwf.dat -c 4 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_cube.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_sum.pov -i ../Analysis/dwf.dat -c 5 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_cube.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_pr.pov   -i ../Analysis/dwf.dat -c 6 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_cube.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_or.pov   -i ../Analysis/dwf.dat -c 7
	# Excecute the python script to generate the dwf profiles along main directions
	$lcDir/dwf_profile.py -p ../w5cb_cube.psf -t $last.dcd -s cube -o profile_dwf_z   -c 4
	$lcDir/dwf_profile.py -p ../w5cb_cube.psf -t $last.dcd -s cube -o profile_dwf_sum -c 5
	$lcDir/dwf_profile.py -p ../w5cb_cube.psf -t $last.dcd -s cube -o profile_dwf_pr   -c 6
	$lcDir/dwf_profile.py -p ../w5cb_cube.psf -t $last.dcd -s cube -o profile_dwf_or   -c 7
	#Navigate to the Analysis directory
	cd $lcDir/cube/Analysis
	#Plot the dwf data using gnuplot
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
	#Prepare the snapshots using POV-ray
	povray +w3000 +H3000 -Iout_dwf_z.pov    +Osnap_dwf_z.png   +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_sum.pov  +Osnap_dwf_sum.png +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_pr.pov    +Osnap_dwf_pr.png   +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_or.pov    +Osnap_dwf_or.png   +Q9 +A0.01 +R9 -D -V
fi
#
####################################################################################################
#squircle @ 25c
if false
then
	first=9
	last=10
	b=216
	e=100000
	#Navigate to the directory
	cd $lcDir/squircle/DCD
	#Excecute the python script to calculate DWF
	#$lcDir/dwf.py -p ../w5cb_squircle.psf -f $first -l $last -nx 100 -ny 100 -b $b -e $e -s squircle
	#Execute the python script to prepare the POV-ray input files
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_squircle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_z.pov   -i ../Analysis/dwf.dat -c 4 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_squircle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_sum.pov -i ../Analysis/dwf.dat -c 5 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_squircle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_pr.pov   -i ../Analysis/dwf.dat -c 6 
	$lcDir/inter_cube.py  -t $last.dcd -p ../w5cb_squircle.psf -l 6 -u 10 -r head_cube -o ../Analysis/out_dwf_or.pov   -i ../Analysis/dwf.dat -c 7
	# Excecute the python script to generate the dwf profiles along main directions
	$lcDir/dwf_profile.py -p ../w5cb_squircle.psf -t $last.dcd -s squircle -o profile_dwf_z   -c 4
	$lcDir/dwf_profile.py -p ../w5cb_squircle.psf -t $last.dcd -s squircle -o profile_dwf_sum -c 5
	$lcDir/dwf_profile.py -p ../w5cb_squircle.psf -t $last.dcd -s squircle -o profile_dwf_pr   -c 6
	$lcDir/dwf_profile.py -p ../w5cb_squircle.psf -t $last.dcd -s squircle -o profile_dwf_or   -c 7
	#Navigate to the Analysis directory
	cd $lcDir/squircle/Analysis
	#Plot the dwf data using gnuplot
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
	#Prepare the snapshots using POV-ray
	povray +w3000 +H3000 -Iout_dwf_z.pov    +Osnap_dwf_z.png   +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_sum.pov  +Osnap_dwf_sum.png +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_pr.pov    +Osnap_dwf_pr.png   +Q9 +A0.01 +R9 -D -V
	povray +w3000 +H3000 -Iout_dwf_or.pov    +Osnap_dwf_or.png   +Q9 +A0.01 +R9 -D -V
fi
#
###################################################################################################

#! /bin/bash

lcDir="/home/ramezani/midway/lc/CUBE"

#
###################################################################################################
#circle @ 25c
if true
then
	first=17
	last=25
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/circle/DCD
	#Excecute the python script to calculate S
	#$lcDir/Q_vector.py -p ../w5cb_circle.psf -f $first -l $last -nx 31 -ny 31 -vx 15 -vy 15 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	#$lcDir/inter_cube.py     -t $last.dcd -p ../w5cb_circle.psf -l 0.3 -u 0.7 -r head_cube -o ../Analysis/out_S.pov   -i ../Analysis/Output_scaler
	# Excecute the python script to generate the density plots
	$lcDir/densityProfile_defectArea.py -p ../w5cb_circle.psf -t $last.dcd -s circle
	#Navigate to "Analysis" directory
	cd $lcDir/circle/Analysis
	#Plot S data using gnuplot
	#gnuplot cc.gnu
	#gnuplot director.gnu
	#Generate image using PovRay
	#povray +w3000 +H3000 -Iout_S.pov    +Osnap_S.png   +Q9 +A0.01 +R9 -D -V
fi	
#
###################################################################################################
#cube @ 25c
if true
then
	first=3
	last=3
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/cube/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	#$lcDir/Q_vector.py -p ../w5cb_cube.psf -f $first -l $last -nx 31 -ny 31 -vx 15 -vy 15 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	#$lcDir/inter_cube.py     -t $last.dcd -p ../w5cb_cube.psf -l 0.3 -u 0.7 -r head_cube -o ../Analysis/out_S.pov   -i ../Analysis/Output_scaler
	# Excecute the python script to generate the density plots
	$lcDir/densityProfile_defectArea.py -p ../w5cb_cube.psf -t $last.dcd -s cube
	#Navigate to "Analysis" directory
	cd $lcDir/cube/Analysis
	#Plot S data using gnuplot
	#gnuplot cc.gnu
	#gnuplot director.gnu
	#Generate image using PovRay
	#povray +w3000 +H3000 -Iout_S.pov    +Osnap_S.png   +Q9 +A0.01 +R9 -D -V
fi
#
####################################################################################################
#squircle @ 25c
if true
then
	first=3
	last=10
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/squircle/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	#$lcDir/Q_vector.py -p ../w5cb_squircle.psf -f $first -l $last -nx 31 -ny 31 -vx 15 -vy 15 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	#$lcDir/inter_cube.py     -t $last.dcd -p ../w5cb_squircle.psf -l 0.3 -u 0.7 -r head_cube -o ../Analysis/out_S.pov   -i ../Analysis/Output_scaler
	# Excecute the python script to generate the density plots
	$lcDir/densityProfile_defectArea.py -p ../w5cb_squircle.psf -t $last.dcd -s squircle
	#Navigate to "Analysis" directory
	cd $lcDir/squircle/Analysis
	#Plot S data using gnuplot
	#gnuplot cc.gnu
	#gnuplot director.gnu
	#Generate image using PovRay
	#povray +w3000 +H3000 -Iout_S.pov    +Osnap_S.png   +Q9 +A0.01 +R9 -D -V
fi
#
###################################################################################################

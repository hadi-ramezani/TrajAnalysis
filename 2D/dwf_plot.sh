#! /bin/bash

lcDir="/home/ramezani/midway/lc/CUBE"

#
###################################################################################################
#circle @ 25c
if true
then
	#Navigate to the Analysis directory
	cd $lcDir/circle/Analysis
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
fi	
#
###################################################################################################
#cube @ 25c
if true
then
	#Navigate to the Analysis directory
	cd $lcDir/cube/Analysis
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
fi
#
####################################################################################################
#squircle @ 25c
if true
then
	#Navigate to the Analysis directory
	cd $lcDir/squircle/Analysis
	gnuplot dwf_radial.gnu
	gnuplot dwf_z.gnu
	gnuplot dwf_sum.gnu
fi
#
###################################################################################################

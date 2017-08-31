#! /bin/bash

lcDir="/home/ramezani/midway/lc/CUBE"

##
####################################################################################################
##circle @ 25c
#first=8
#last=25
#b=0
#e=100000
##Navigate to the directory
#cd $lcDir/circle/DCD
##Excecute the python script to calculate S
#$lcDir/anchoring.py  -p ../w5cb_circle.psf -f $first -l $last -n 320 -T 294.15 -s circle -b $b -e $e 
##
###################################################################################################
#cube @ 25c
first=2
last=3
b=0
e=100000
#Navigate to the directory
cd $lcDir/cube/DCD
#Excecute the python script to calculate S and P2 and cos(beta) Profile
$lcDir/anchoring.py  -p ../w5cb_cube.psf -f $first -l $last -n 320 -T 294.15 -s cube -b $b -e $e 
#
###################################################################################################
#squircle @ 25c
first=2
last=10
b=0
e=100000
#Navigate to the directory
cd $lcDir/squircle/DCD
#Excecute the python script to calculate S and P2 and cos(beta) Profile
$lcDir/anchoring.py  -p ../w5cb_squircle.psf -f $first -l $last -n 320 -T 294.15 -s squircle -b $b -e $e 
#
###################################################################################################

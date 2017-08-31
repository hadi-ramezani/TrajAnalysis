#! /bin/bash

#lcDir="/home/ramezani/midway/lc"
lcDir="/project/roux/ramezani/lc"
#
###################################################################################################
#5CB_vac @ 16c
if false
then
	first=73
	last=80
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/vac/16c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.5 -u 0.8 -r head_5cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.3 -u 0.3 -r head_5cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l  0.0 -u 1.0 -r head_5cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/vac/16c/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov    +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov  +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov    +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_vac @ 20c
if false
then 
	first=73
	last=80
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/vac/20c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.5 -u 0.8 -r head_5cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat 
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.3 -u 0.3 -r head_5cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat 
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l  0.0 -u 1.0 -r head_5cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/vac/20c/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_vac @ 25c
if false
then
	first=73
	last=80
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/vac/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_vac.psf -f $first -l $last 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.5 -u 0.8 -r head_5cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.3 -u 0.3 -r head_5cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l  0.0 -u 1.0 -r head_5cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/vac/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_vac @ 38c
if false
then
	first=73
	last=80
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/vac/38c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.5 -u 0.8 -r head_5cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l -0.3 -u 0.3 -r head_5cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w5cb_vac.psf -l  0.0 -u 1.0 -r head_5cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/vac/38c/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_vac @ 32c
if false
then
	first=6
	last=6
	b=5500
	e=7000
	#Navigate to the directory
	cd $lcDir/8cb/vac/32c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.5 -u 0.8 -r head_8cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.3 -u 0.3 -r head_8cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l  0.0 -u 1.0 -r head_8cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/vac/32c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_vac @ 38c, Note: first=83, last=83, b=7000, e=7499
if false
then
	first=87
	last=87
	b=2500
	e=2900
	#Navigate to the directory
	cd $lcDir/8cb/vac/38c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.5 -u 0.8 -r head_8cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.3 -u 0.3 -r head_8cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l  0.0 -u 1.0 -r head_8cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/vac/38c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_vac @ 43c
if false
then
	first=73
	last=80
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/vac/43c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.5 -u 0.8 -r head_8cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l -0.3 -u 0.3 -r head_8cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w8cb_vac.psf -l  0.0 -u 1.0 -r head_8cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/vac/43c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_wat @ 25c
if false
then
	first=77
	last=78
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/wat/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/wat/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_wat @ 40c
if false
then
	first=90
	last=97
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/wat/40c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w5cb_wat.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/wat/40c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
##################################################################################################
#5CB_wat_nacl @ 25c
if false
then
	first=71
	last=71
	b=8000
	e=9499
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name NY1 and prop z < 320" -o ../Analysis/dist_N.dat      -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name OH2"                  -o ../Analysis/dist_O.dat      -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name SOD"                  -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name CLA"                  -o ../Analysis/dist_anion.dat  -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nacl.py
	Navigate to "Analysis" directory
	cd $lcDir/5cb/nacl/2m/25c/Analysis
	Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_wat_nacl @ 40c
if false
then
	first=97
	last=104
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/40c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w5cb_wat_nacl.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 1 -s "name CLA" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nacl.py
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/nacl/2m/40c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_wat_nai @ 25c
if false
then
	first=78
	last=78
	b=4500
	e=5999
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name NY1"                  -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name OH2"                  -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name SOD"                  -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name I"                  -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nai.py
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/nai/2m/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#5CB_wat_nai @ 40c
if false
then
	first=94
	last=101
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/40c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l -0.5 -u 0.8 -r head_5cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l -0.3 -u 0.3 -r head_5cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w5cb_wat_nai.psf -l  0.0 -u 1.0 -r head_5cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 1 -s "name I" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nai.py
	#Navigate to "Analysis" directory
	cd $lcDir/5cb/nai/2m/40c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3000 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3000 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat @ 25c
if false
then
	first=1
	last=1
	b=8750
	e=9000
	#Navigate to the directory
	cd $lcDir/8cb/wat/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/wat/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_32 @ 32c 
if false
then
	first=4
	last=4
	b=7800
	e=8344
	#Navigate to the directory
	cd $lcDir/8cb/wat/32c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/wat/32c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat @ 37c Note: first=8, last=8, b=9750, e=100000
# Note first=12, last=12, b=7500, e=8019
if false
then
	first=15
	last=15
	b=7850
	e=8391
	#Navigate to the directory
	cd $lcDir/8cb/wat/37c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/wat/37c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#
#8CB_wat_big @ 37c 
if true
then
	first=28
	last=28
	b=0
	e=100000
	#Navigate to the directory
	cd /project/depablo/ramezani/lc/8cb/wat/big_37c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 200 -sc "name CA12 and prop z < 860" -sn "name NY1 and prop z < 860" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 200 -sc "name CA12 and prop z < 860" -sn "name NY1 and prop z < 860" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd /project/depablo/ramezani/lc/8cb/wat/big_37c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat @ 47c
if false
then
	first=70
	last=74
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/wat/47c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_wat.py     -t $last.dcd -p ../w8cb_wat.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_wat.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/wat/47c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nacl @ 25c
if false
then
	first=44
	last=46
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name CLA" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nacl.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nacl/2m/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nacl @ 37c
if false
then
	first=44
	last=46
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/37c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name CLA" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nacl.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nacl/2m/37c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nacl @ 47c
if false
then
	first=78
	last=82
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/47c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nacl.py     -t $last.dcd -p ../w8cb_wat_nacl.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 1 -s "name CLA" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nacl.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nacl/2m/47c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nai @ 25c
if false
then
	first=39
	last=41
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/25c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name I" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nai.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nai/2m/25c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nai @ 37c
if false
then
	first=39
	last=41
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/37c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name I" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nai.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nai/2m/37c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#8CB_wat_nai @ 47c
if false
then
	first=72
	last=76
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/47c/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	$lcDir/tensorZ.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.5 -u 0.8 -r head_8cb_wat -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l -0.3 -u 0.3 -r head_8cb_wat -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_nai.py     -t $last.dcd -p ../w8cb_wat_nai.psf -l  0.0 -u 1.0 -r head_8cb_wat -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Execute the python script to calculate distribution of atoms
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name NY1" -o ../Analysis/dist_N.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name OH2" -o ../Analysis/dist_O.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name SOD" -o ../Analysis/dist_cation.dat -b $b -e $e
	$lcDir/dist.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 1 -s "name I" -o ../Analysis/dist_anion.dat -b $b -e $e
	#Execute the python script to plot the distribution data
	$lcDir/dist_plot_nai.py
	#Navigate to "Analysis" directory
	cd $lcDir/8cb/nai/2m/47c/Analysis
	#Generate image using PovRay
	povray +w1000 +H3500 -Iout_p.pov   +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_cos.pov +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H3500 -Iout_s.pov   +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#50_8cb_50_5cb @ 38c
if false
then
	first=8
	last=11
	b=2000
	e=100000
	#Navigate to the directory
	cd $lcDir/mix/50_8cb_50_5cb/38/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l -0.5 -u 0.8 -r head_mix_50 -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l -0.3 -u 0.3 -r head_mix_50 -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l  0.0 -u 1.0 -r head_mix_50 -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/mix/50_8cb_50_5cb/38/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov    +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov  +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov    +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#50_8cb_50_5cb @ 27c
if false
then
	first=3
	last=3
	b=16000
	e=100000
	#Navigate to the directory
	cd $lcDir/mix/50_8cb_50_5cb/27/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l -0.5 -u 0.8 -r head_5cb -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l -0.3 -u 0.3 -r head_5cb -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w50_8cb_50_5cb.psf -l  0.0 -u 1.0 -r head_5cb -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/mix/50_8cb_50_5cb/27/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov    +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov  +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov    +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#85_8cb_15_5cb @ 38c Note: previous settings: first=13, last=15, b=1700, e=100000
if false
then
	first=16
	last=16
	b=14200
	e=100000
	#Navigate to the directory
	cd $lcDir/mix/85_8cb_15_5cb/38/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l -0.5 -u 0.8 -r head_mix_85 -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l -0.3 -u 0.3 -r head_mix_85 -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l  0.0 -u 1.0 -r head_mix_85 -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/mix/85_8cb_15_5cb/38/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov    +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov  +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov    +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#
###################################################################################################
#85_8cb_15_5cb @ 27c
if false
then
	first=3
	last=3
	b=14200
	e=100000
	#Navigate to the directory
	cd $lcDir/mix/85_8cb_15_5cb/27/DCD
	#Excecute the python script to calculate S and P2 and cos(beta) Profile
	$lcDir/sZ_cosZ.py  -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	$lcDir/tensorZ.py  -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e
	#Execute the python script to generate the input file for PovRay
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l -0.5 -u 0.8 -r head_mix_85 -o ../Analysis/out_p.pov   -i ../Analysis/ScalerZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l -0.3 -u 0.3 -r head_mix_85 -o ../Analysis/out_cos.pov -i ../Analysis/cosZ.dat
	$lcDir/inter_vac.py     -t $last.dcd -p ../w85_8cb_15_5cb.psf -l  0.0 -u 1.0 -r head_mix_85 -o ../Analysis/out_s.pov   -i ../Analysis/tensorZ.dat_local
	#Navigate to "Analysis" directory
	cd $lcDir/mix/85_8cb_15_5cb/27/Analysis
	#Generate image using PovRay
	povray +w1000 +H2500 -Iout_p.pov    +Osnap_p.png   +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_cos.pov  +Osnap_cos.png +Q10 +A0.01 +R9 -D -V
	povray +w1000 +H2500 -Iout_s.pov    +Osnap_s.png   +Q10 +A0.01 +R9 -D -V
fi
#

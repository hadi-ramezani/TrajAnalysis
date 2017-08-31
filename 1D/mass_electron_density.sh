#! /bin/bash

lcDir="/home/ramezani/midway/lc"

#
###################################################################################################
#5CB_vac @ 16c
if false
then
	first=103
	last=103
	b=0
	e=5000
	#Navigate to the directory
	cd $lcDir/5cb/vac/16c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#5CB_vac @ 20c
if false
then
	first=104
	last=104
	b=0
	e=5000
	#Navigate to the directory
	cd $lcDir/5cb/vac/20c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#5CB_vac @ 25c
if false
then
	first=104
	last=104
	b=0
	e=5000
	#Navigate to the directory
	cd $lcDir/5cb/vac/25c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#5CB_vac @ 38c
if false
then
	first=105
	last=105
	b=0
	e=5000
	#Navigate to the directory
	cd $lcDir/5cb/vac/38c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_vac.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#8CB_vac @ 38c
if false
then
	first=83
	last=83
	b=7000
	e=7499
	#Navigate to the directory
	cd $lcDir/8cb/vac/38c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_vac.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_vac.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
####################################################################################################
#8CB_vac @ 43c
if false
then
	first=123
	last=123
	b=0
	e=5000
	#Navigate to the directory
	cd $lcDir/8cb/vac/43c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_vac.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_vac.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w5cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 5CB and prop z < 320" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#8CB_wat @ 25c
if false
then
	first=36
	last=38
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/wat/25c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################
#8CB_wat @ 37c
if false
then
	first=35
	last=37
	b=0
	e=100000
	#Navigate to the directory
	cd $lcDir/8cb/wat/37c/DCD
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
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
	#Execute the python script to calculate mass density of LC phase and plot the data
	$lcDir/mass_density.py     -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	$lcDir/electron_density.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 3 -s "resname 8CB and prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
fi
#
###################################################################################################                                                  
#50_8cb_50_5cb @ 38c                                                                                                                                 
if true                                                                                                                                             
then                                                                                                                                                 
	first=8                                                                                                                                             
	last=11                                                                                                                                              
	b=2000                                                                                                                                           
	e=100000                                                                                                                                         
	#Navigate to the directory                                                                                                                          
	cd $lcDir/mix/50_8cb_50_5cb/38/DCD                                                                                                                  
	#Execute the python script to calculate mass density of LC phase and plot the data
	#$lcDir/mass_density.py     -p ../w50_8cb_50_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	#$lcDir/electron_density.py -p ../w50_8cb_50_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
	$lcDir/mix_density.py -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e -p5 50 -p8 50
fi                                                                                                                                                   
#                                                                                                                                                    
###################################################################################################                                                  
#50_8cb_50_5cb @ 27c                                                                                                                                 
if true                                                                                                                                            
then                                                                                                                                                 
	first=3                                                                                                                                             
	last=3                                                                                                                                              
	b=16000                                                                                                                                             
	e=100000                                                                                                                                            
	#Navigate to the directory                                                                                                                          
	cd $lcDir/mix/50_8cb_50_5cb/27/DCD                                                                                                                  
	#Execute the python script to calculate mass density of LC phase and plot the data
	#$lcDir/mass_density.py     -p ../w50_8cb_50_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	#$lcDir/electron_density.py -p ../w50_8cb_50_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
	$lcDir/mix_density.py -p ../w50_8cb_50_5cb.psf -f $first -l $last -n 100 -b $b -e $e -p5 50 -p8 50
fi                                                                                                                                                   
#                                                                                                                                                    
###################################################################################################                                                  
#85_8cb_15_5cb @ 38c                                                                                                                                 
if true                                                                                                                                             
then                                                                                                                                                 
	first=13                                                                                                                                            
	last=15                                                                                                                                             
	b=1700                                                                                                                                                 
	e=100000                                                                                                                                            
	#Navigate to the directory                                                                                                                          
	cd $lcDir/mix/85_8cb_15_5cb/38/DCD                                                                                                                  
	#Execute the python script to calculate mass density of LC phase and plot the data
	#$lcDir/mass_density.py     -p ../w85_8cb_15_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	#$lcDir/electron_density.py -p ../w85_8cb_15_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
	$lcDir/mix_density.py -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e -p5 15 -p8 85
fi                                                                                                                                                   
#                                                                                                                                                    
###################################################################################################                                                  
#85_8cb_15_5cb @ 27c                                                                                                                                 
if true                                                                                                                                             
then                                                                                                                                                 
	first=3                                                                                                                                             
	last=3                                                                                                                                             
	b=14200                                                                                                                                                 
	e=100000                                                                                                                                            
	#Navigate to the directory                                                                                                                          
	cd $lcDir/mix/85_8cb_15_5cb/27/DCD                                                                                                                  
	#Execute the python script to calculate mass density of LC phase and plot the data
	#$lcDir/mass_density.py     -p ../w85_8cb_15_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/mass_density.dat      -b $b -e $e
	#$lcDir/electron_density.py -p ../w85_8cb_15_5cb.psf -f $first -l $last -d 3 -s "prop z < 450" -o ../Analysis/electron_density.dat  -b $b -e $e
	$lcDir/mix_density.py -p ../w85_8cb_15_5cb.psf -f $first -l $last -n 100 -b $b -e $e -p5 15 -p8 85
fi                                                                                                                                                   
#                                                                                                                                                    
                                                                                                                                                     
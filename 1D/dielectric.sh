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
	#Navigate to the directory
	cd $lcDir/5cb/vac/16c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w5cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 20c
if false
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/20c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w5cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 25c
if false
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/25c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w5cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 38c
if false
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/38c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w5cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#8CB_vac @ 38c
if false
then
	first=87
	last=87
	#Navigate to the directory
	cd $lcDir/8cb/vac/38c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w8cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#8CB_vac @ 43c
if false
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/8cb/vac/43c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/dielectric.py  -p ../w8cb_vac.psf -f $first -l $last -d 5 -s "name CA12" -o ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat @ 25c
if false
then
	first=70
	last=78
	#Navigate to the directory
	cd $lcDir/5cb/wat/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat.psf -f $first -l $last -d 2 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat.psf -f $first -l $last -d 2 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat @ 40c
if false
then
	first=70
	last=97
	#Navigate to the directory
	cd $lcDir/5cb/wat/40c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat.psf -f $first -l $last -d 5 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat.psf -f $first -l $last -d 5 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
##################################################################################################
#5CB_wat_nacl @ 25c
if true
then
	first=70
	last=71
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 5 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 5 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nacl @ 40c
if true
then
	first=70
	last=104
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/40c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 5 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 5 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nai @ 25c
if true
then
	first=70
	last=78
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 5 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 5 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nai @ 40c
if true
then
	first=70
	last=101
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/40c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/dielectric.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 5 -s "name CA12 and prop z < 320" -o ../Analysis/lc_dielectric.dat                                   #-b 1 -e 5
	$lcDir/dielectric.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 5 -s "name OH2"                   -o ../Analysis/wat_dielectric.dat -i ../Analysis/lc_dielectric.dat #-b 1 -e 5
fi
#
###################################################################################################
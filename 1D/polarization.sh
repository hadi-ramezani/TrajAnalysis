#! /bin/bash

lcDir="/home/ramezani/midway/lc"

#
###################################################################################################
#5CB_vac @ 16c
if true
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/16c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w5cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 20c
if true
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/20c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w5cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 25c
if true
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/25c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w5cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_vac @ 38c
if true
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/5cb/vac/38c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w5cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_vac @ 38c
if true
then
	first=87
	last=87
	#Navigate to the directory
	cd $lcDir/8cb/vac/38c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w8cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_vac @ 43c
if true
then
	first=73
	last=80
	#Navigate to the directory
	cd $lcDir/8cb/vac/43c/DCD
	#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
	$lcDir/polarization_lc_only.py  -p ../w8cb_vac.psf -f $first -l $last -d 2 #-b 1 -e 5
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
	$lcDir/polarization.py -p ../w5cb_wat.psf -f $first -l $last -d 2 #-b 1 -e 5
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
	$lcDir/polarization.py -p ../w5cb_wat.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
##################################################################################################
#5CB_wat_nacl @ 25c
if false
then
	first=70
	last=71
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nacl @ 40c
if false
then
	first=70
	last=104
	#Navigate to the directory
	cd $lcDir/5cb/nacl/2m/40c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w5cb_wat_nacl.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nai @ 25c
if false
then
	first=70
	last=78
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#5CB_wat_nai @ 40c
if false
then
	first=70
	last=101
	#Navigate to the directory
	cd $lcDir/5cb/nai/2m/40c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w5cb_wat_nai.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat @ 25c
if false
then
	first=30
	last=38
	#Navigate to the directory
	cd $lcDir/8cb/wat/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat @ 37c
if false
then
	first=30
	last=37
	#Navigate to the directory
	cd $lcDir/8cb/wat/37c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat @ 47c
if false
then
	first=30
	last=74
	#Navigate to the directory
	cd $lcDir/8cb/wat/47c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nacl @ 25c
if false
then
	first=30
	last=46
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nacl @ 37c
if false
then
	first=30
	last=46
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/37c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nacl @ 47c
if false
then
	first=30
	last=82
	#Navigate to the directory
	cd $lcDir/8cb/nacl/2m/47c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nacl.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nai @ 25c
if false
then
	first=30
	last=41
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/25c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nai @ 37c
if false
then
	first=30
	last=41
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/37c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
###################################################################################################
#8CB_wat_nai @ 47c
if false
then
	first=30
	last=76
	#Navigate to the directory
	cd $lcDir/8cb/nai/2m/47c/DCD
	#Execute the python script to calculate distribution of atoms
	$lcDir/polarization.py -p ../w8cb_wat_nai.psf -f $first -l $last -d 2 #-b 1 -e 5
fi
#
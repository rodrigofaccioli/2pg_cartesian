#!/bin/sh
# This script runs g_energy program only.

# $1 where Gromacs programs are.
# $2 where local execute is. Here is your working directory where all files will be created
# $3 The value to get the energy. 
# $4 File name where the energy is stored.

#g_energy - prepare the file to read the potential energy
# echo $3 | $1./g_energy -f $2file_energy_computed.ener.edr -o $2$4 > /dev/null 2> /dev/null

if [ "$1" = "g_energy" ]; then
	echo "$5" | "$2"./g_energy -f "$3" -o "$4" #> /dev/null 2> /dev/null
elif [ "$1" = 'g_sas' ]; then
	# g_sas - compute areas
	echo "$6" echo "$7" | "$2"./g_sas -f "$3" -s "$4" -o "$5" #> /dev/null 2> /dev/null
elif [ "$1" = 'g_gyrate' ]; then
	echo "$6" | "$2"./g_gyrate -f "$3" -s "$4" -o "$5" #> /dev/null 2> /dev/null
elif [ "$1" = "g_hbond" ]; then
	echo "$6" echo "$7" | "$2"./g_hbond -f "$3" -s "$4" -num "$5" #> /dev/null 2> /dev/null
fi


#echo "$4" | "$1"./g_energy -f "$2" -o "$3"

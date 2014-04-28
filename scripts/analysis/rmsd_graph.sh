#!/bin/bash
# This script makes a plot of the best RMSD by generation
#
# $1 where Gromacs programs are.
# $2 where local execute is. Here is your working directory where all files will be created
# $3 protpred executables location.
# $4 Configuration file
# $5 Number of generations
# $6 file in which the best RMSD AlfaCarbon-AlfaCarbon of each generation will be stored
# $7 file in which the best RMSD of backbone-backbone each generation will be stored
# $8 file in which the best RMSD of AlfaCarbon-Side-Chain each generation will be stored
# $9 file in which the best RMSD of Side-Chain-Side-Chain each generation will be stored

# counter variable
i=1

# If the output files exists, delete them
if [ -e $2$6 ]; then rm $2$6; fi
if [ -e $2$7 ]; then rm $2$7; fi
if [ -e $2$8 ]; then rm $2$8; fi
if [ -e $2$9 ]; then rm $2$9; fi
if [ -e $2temp_RMSD.rmsd ]; then rm $2temp_RMSD.rmsd; fi

# Creates empty output files
touch $2$6 $2$7

let total=$5+1

while [ $i -lt $total ]; do

        # Generates pdb files for the ith generation
        $3./protpred-Gromacs_Nerf_Population $2$4 pop_"$i".pop

        # Calculates and classifies individuals according to their RMSD
        $3./protpred-Gromacs-RMSD $2$4

        # Adds the best RMSD to the output file
        head -n 2 protpred-Gromacs.rmsd | tail -n 1 >> $2$6
        head -n 2 protpred-Gromacs_backbone.rmsd | tail -n 1 >> $2$7
        head -n 2 protpred-Gromacs_alfaCarbon-SideChain.rmsd | tail -n 1 >> $2$8
        head -n 2 protpred-Gromacs_SideChain-SideChain.rmsd | tail -n 1 >> $2$9

	# Adds current generation to the checking file. useful when reading with "tail -f"
        echo $i >> $2temp_RMSD.rmsd

	# Increments the counter
        let i=$i+1

done

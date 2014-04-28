#!/bin/bash
# This script makes a pdb trajectory containing the best individual of each generation
#
# $1 where Gromacs programs are.
# $2 where local execute is. Here is your working directory where all files will be created
# $3 protpred executables location.
# $4 Configuration file
# $5 Number of generations
# $6 file in which a trajectory containing the lowest RMSD individual for each generation will be created
# $7 file for checking the script current generation


# Counter variable
i=1

# If the output files exists, delete them
if [ -e $2$6 ]; then rm $2$6; fi
if [ -e $2$7 ]; then rm $2$7; fi

# Creates empty output files
touch $2$6 $2$7

let total=$5+1

while [ $i -lt $total ]; do

        # Generates pdb files for the ith generation
        $3./protpred-Gromacs_Nerf_Population $2$4 pop_"$i".pop

        # Calculates and classifies individuals according to their RMSD
        $3./protpred-Gromacs-RMSD $2$4

	# spaces.rmsd is a necessary file because of the spaces between the columns of the original file
        more protpred-Gromacs.rmsd | sed 's/\ /\ \ \ \ \ /g' > spaces.rmsd

        # Identifies the best individues
        best=$(head -n 5 spaces.rmsd | tail -n 1 | cut -c 1-22 | sed 's/\ //g')

        # Copies the best individual pdb
        cp $best rmsd_temporary_"$i".pdb

	# Adds the MODEL tag to the output pdb file
        echo MODEL $i >> $2$6

	# Adds the coordinates to the output file
        more rmsd_temporary_"$i".pdb | grep ATOM >> $2$6

	# Adds the ENDMDL tag to the output file
        echo ENDMDL >> $2$6

	# Adds current generation to the checking file. useful when reading with "tail -f"
        echo $i >> $2$7

        let i=$i+1

done

# Adds the END tag to the output file
echo END >> $2$6

# Fits the trajectory
echo 3 1 | trjconv -f $2$6 -o $2$6 -fit progressive -s $2$6

rm spaces.rmsd rmsd_temporary_*

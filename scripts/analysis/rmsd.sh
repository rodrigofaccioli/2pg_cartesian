#!/bin/bash

#
#
#	Script para calcular RMSD dos indivíduos não dominados.
#
#							Leandro Oliveira Bortot
#							26 Maio 2014
#							Rodrigo Antonio Faccioli
#							27 Maio 2014
#
#

total_gen=$1
native=$2
path_2pg=$3
path_gmx=$4


# Apaga arquivos pre-existentes
rm "plot_non-dominated_rmsd_avg-sd.xvg" 2>/dev/null
rm "plot_non-dominated_rmsd_scatter.xvg" 2>/dev/null
rm "plot_non-dominated_rmsd_minimum.xvg" 2>/dev/null
rm "temporary_"* 2>/dev/null



gen=1
while [ $gen -le $total_gen ]; do	# For each generation ... 

# RMSD
	# Calculates the Ca-Ca RMSD of the non-dominated individuals
	echo "C-alpha C-alpha" | "$path_gmx"./g_rms -f "pop_NON_DOMINATED_""$gen"".pdb" -s "$native" -o "temporary_rmsd.xvg"
	
	# Take the rmsd values
	grep -v "#" "temporary_rmsd.xvg" | grep -v "@" | awk '{print $2}' | sort -g > "temporary_rmsd_values"
	
	# Calculates the average and standard deviation
	avgsd=$("$path_2pg""scripts/analysis/"./avg-sd "temporary_rmsd_values" $total_ind_front)
	
	# Plot
	echo "$gen"" ""$avgsd" >> "plot_non-dominated_rmsd_avg-sd.xvg"
	
	
	
	# Plots the smaller value for the current generation
	head -n 1 "temporary_rmsd_values" | awk '{print '"$gen"',$1}' >> "plot_non-dominated_rmsd_minimum.xvg"
	
	
	
	# Reads the rmsd values and make the scatter plot
	grep -v "#" "temporary_rmsd.xvg" | grep -v "@" | awk '{print '"$gen"',$2}' >> "plot_non-dominated_rmsd_scatter.xvg"
	
	rm \#* 2>/dev/null

	let gen=$gen+1
	
done


rm "temporary_"*
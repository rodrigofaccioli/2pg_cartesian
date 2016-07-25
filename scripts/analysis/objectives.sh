#!/bin/bash

#
#
#	Script para analisar predições com algoritmo evolutivo. São analisados:
#		- O valor de cada objetivo dos indivíduos não-dominados ao longo das gerações
#		- Número de indivíduos não-dominados por geração
#
#							Leandro Oliveira Bortot
#							26 Maio 2014
#
#

total_gen=$1
path_2pg=$2
path_gmx=$3

# Apaga arquivos pre-existentes
rm "plot_non-dominated_"*"_avg-sd.xvg" 2>/dev/null
rm "plot_non-dominated_"*"_scatter.xvg" 2>/dev/null
rm "plot_non-dominated_nfront.xvg" 2>/dev/null
rm "temporary_"* 2>/dev/null

# Creates a list with the names of all objectives
find . -maxdepth 1 -type f -name "*_NON_DOMINATED_1.fit" | sed 's/.\///g' | sed 's/_NON_DOMINATED_1.fit//g' > "temporary_obj_list"

# Total number of objectives
total_obj=$(wc -l "temporary_obj_list" | awk '{print $1}')



gen=1
while [ $gen -le $total_gen ]; do	# For each generation ... 

	echo -e "\tGeneration $gen"" / ""$total_gen"

	
	obj=1
	while [ $obj -le $total_obj ]; do	# For each objective ...
	
		obj_name=$(head -n $obj "temporary_obj_list" | tail -n 1 | awk '{print $1}')	# Gets the name of the current objective
		
		
		
# PLOTS THE OBJECTIVE VALUES BY GENERATION
		# Reads the .fit file for the current objective and generation and get the objective values
		grep -v "#" "$obj_name""_NON_DOMINATED_""$gen"".fit" | grep -v "@" | awk '{print $2}' > "temporary_obj_values"
		
		# Total of non-dominated individuals in this generation
		total_ind_front=$(wc -l "temporary_obj_values" | awk '{print $1}')
		
		# Calculates the average and standard deviation of the values
		avgsd=$("$path_2pg""scripts/analysis/"./avg-sd "temporary_obj_values" $total_ind_front)
		
		# Plot
		echo "$gen"" ""$avgsd" >> "plot_non-dominated_""$obj_name""_avg-sd.xvg"
		
		
		
		# Reads the .fit files and plot the objective values for the current generation and objective
		grep -v "#" "$obj_name""_NON_DOMINATED_""$gen"".fit" | grep -v "@" | awk '{print '"$gen"',$2}' >> "plot_non-dominated_""$obj_name""_scatter.xvg"
		
		
		let obj=$obj+1
		
	done
	
		
		
# PLOTS THE NUMBER OF NON-DOMINATED INDIVIDUALS
	grep -v "#" "$obj_name""_NON_DOMINATED_""$gen"".fit" | grep -v "@" | wc -l | awk '{print '"$gen"',$1}' >> "plot_non-dominated_nfront.xvg"
	
	
	let gen=$gen+1
done

	
	
	
#!/bin/bash

find . -maxdepth 1 -type f -name "all_*.front" | sed 's/\.\///g' > "temporary_front_list"

total_front=$(wc -l "temporary_front_list" | awk '{print $1}')

# Makes the individual list
front_name=$(head -n 1 "temporary_front_list" | tail -n 1)
grep -v "#" "$front_name" | awk '{print $6}' > "temporary_ind_list"
total_ind=$(wc -l "temporary_ind_list" | awk '{print $1}')


echo -e "avg-ranking""\t""stdev-ranking""\t""avg-front""\t""stdev-front""\t""avg-dominated""\t""stdev-dominated""\t""method" > "final.dat"

ind=1
while [ $ind -le $total_ind ]; do
	
	ind_name=$(head -n $ind "temporary_ind_list" | tail -n 1)
	
	echo "$ind_name"
#	front=1
#	while [ $front -le $total_front ]; do
	
#		front_name=$(head -n $front "temporary_front_list" | tail -n 1)
	
		grep "$ind_name" "all_front_"*".front" | awk '{print $2}' > "temporary_ranking"
		grep "$ind_name" "all_front_"*".front" | awk '{print $3}' > "temporary_front"
		grep "$ind_name" "all_front_"*".front" | awk '{print $4}' > "temporary_dominated"
		
		~/d/scripts/analise/avg-sd "temporary_ranking" $total_front > "temporary_average_ranking"
		~/d/scripts/analise/avg-sd "temporary_front" $total_front > "temporary_average_front"
		~/d/scripts/analise/avg-sd "temporary_dominated" $total_front > "temporary_average_dominated"
		
		echo "$ind_name" > "temporary_name"
		
		paste "temporary_average_ranking" "temporary_average_front" "temporary_average_dominated" "temporary_name" >> "temporary_results"
		
		
		#sed 's/\ //g' "$front_name" | grep -v $'\t'"$front"$'\t\t' | grep $'\t'"$front"$'\t'
	
		#let front=$front+1
	#done
		
	let ind=$ind+1
done


sort -g "temporary_results" >> "final.dat"


rm "temporary_"*
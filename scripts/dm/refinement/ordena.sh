#!/bin/bash

search=$1

rm "temporary_xvg_list" "temporary_values" "temporary_names" "temporary_list" 2>/dev/null



find . -type f -name "*""$search""*.xvg" | sed 's/.\///g' >> "temporary_xvg_list"


total=$(wc -l "temporary_xvg_list" | awk '{print $1}')

i=1
while [ $i -le $total ]; do

	name=$(head -n $i "temporary_xvg_list" | tail -n 1)

	tail -n 1 "$name" | awk '{print $2}' >> "temporary_values"
	echo "$name" >> "temporary_names"

	let i=$i+1

done


paste "temporary_values" "temporary_names" >> "temporary_list"
sort -g "temporary_list" > "sorted_""$search"".dat"


rm "temporary_xvg_list" "temporary_values" "temporary_names" "temporary_list"

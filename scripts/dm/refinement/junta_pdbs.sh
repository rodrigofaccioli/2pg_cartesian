#!/bin/bash

output=$1

find . -type f | sed 's/.\///g' > "temporary_pdb_list"


total=$(wc -l "temporary_pdb_list" | awk '{print $1}')

i=1
while [ $i -le $total ]; do

	name=$(head -n $i "temporary_pdb_list" | tail -n 1)

	echo "MODEL ""$i" >> "$output"
	cat "$name" | grep "ATOM" >> "$output"
	echo "ENDMDL" >> "$output"

	let i=$i+1

done
rm "temporary_pdb_list"

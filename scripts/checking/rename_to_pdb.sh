#!/bin/bash

# ***** WARNING ALL FILES WILL BE RENAMED  ****************

#path where files are
path=$1

here=$(pwd)

cd "$path"
find . -type f | sed 's/.\///g' | grep -v "list" > "list"

total=$(wc -l "list" | awk '{print $1}')


i=1
while [ $i -le $total ]; do

	name=$(head -n $i "list" | tail -n 1)

	echo "$name"

	mv "$name" "$name"".pdb"

	let i=$i+1
done

rm "list"
cd "$here"

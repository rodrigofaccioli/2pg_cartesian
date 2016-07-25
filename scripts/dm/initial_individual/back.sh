#!/bin/bash

here=$(pwd)
from=$1
to=$2

mkdir -p "$to"
cd "$to"

n=1
while [ $n -le 10 ]; do

	echo -n "vm""$n"" ..."

	ssh vm"$n" "cd ""$from"" ; tar -czf ""$n"".tar.gz md.*"
	scp vm"$n":"$from""/""$n"".tar.gz" .
	tar -xzf "$n"".tar.gz"

	echo "ok"

	let n=$n+1

done



cd "$here"

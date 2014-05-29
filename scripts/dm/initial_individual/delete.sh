#!/bin/bash

directory=$1

n=1
while [ $n -le 10 ]; do

	echo -n "vm""$n"" ..."

	ssh vm"$n" "rm -r ""$directory"

	echo "ok"

	let n=$n+1
done

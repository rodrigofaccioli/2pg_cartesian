#!/bin/bash

nodefile=$1
path_to_pdbs=$2
destination=$3

here=$(pwd)
total_rep=5


total_nodes=$(wc -l "$nodefile" | awk '{print $1}')

node=1
while [ $node -le $total_nodes ]; do

	node_name=$(head -n $node "$nodefile" | tail -n 1)

	echo -n "$node_name"" ... "

	scp "$node_name":"$path_to_pdbs""/"* "$destination" >/dev/null 

	echo "ok"

	let node=$node+1
done
